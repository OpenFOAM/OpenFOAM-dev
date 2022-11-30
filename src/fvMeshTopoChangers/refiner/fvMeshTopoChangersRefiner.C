/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "fvMeshTopoChangersRefiner.H"
#include "surfaceInterpolate.H"
#include "polyTopoChange.H"
#include "syncTools.H"
#include "pointFields.H"
#include "sigFpe.H"
#include "cellSet.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshTopoChangers
{
    defineTypeNameAndDebug(refiner, 0);
    addToRunTimeSelectionTable(fvMeshTopoChanger, refiner, fvMesh);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::label Foam::fvMeshTopoChangers::refiner::count
(
    const PackedBoolList& l,
    const unsigned int val
)
{
    label n = 0;
    forAll(l, i)
    {
        if (l.get(i) == val)
        {
            n++;
        }

        // Debug also serves to get-around Clang compiler trying to optimise
        // out this forAll loop under O3 optimisation
        if (debug)
        {
            Info<< "n=" << n << endl;
        }
    }

    return n;
}


void Foam::fvMeshTopoChangers::refiner::calculateProtectedCells
(
    PackedBoolList& unrefineableCells
) const
{
    if (protectedCells_.empty())
    {
        unrefineableCells.clear();
        return;
    }

    const labelList& cellLevel = meshCutter_.cellLevel();

    unrefineableCells = protectedCells_;

    // Get neighbouring cell level
    labelList neiLevel(mesh().nFaces() - mesh().nInternalFaces());

    for
    (
        label facei = mesh().nInternalFaces();
        facei < mesh().nFaces();
        facei++
    )
    {
        neiLevel[facei - mesh().nInternalFaces()] =
            cellLevel[mesh().faceOwner()[facei]];
    }
    syncTools::swapBoundaryFaceList(mesh(), neiLevel);


    while (true)
    {
        // Pick up faces on border of protected cells
        boolList seedFace(mesh().nFaces(), false);

        forAll(mesh().faceNeighbour(), facei)
        {
            const label own = mesh().faceOwner()[facei];
            const bool ownProtected = unrefineableCells.get(own);
            const label nei = mesh().faceNeighbour()[facei];
            const bool neiProtected = unrefineableCells.get(nei);

            if (ownProtected && (cellLevel[nei] > cellLevel[own]))
            {
                seedFace[facei] = true;
            }
            else if (neiProtected && (cellLevel[own] > cellLevel[nei]))
            {
                seedFace[facei] = true;
            }
        }

        for
        (
            label facei = mesh().nInternalFaces();
            facei < mesh().nFaces();
            facei++
        )
        {
            const label own = mesh().faceOwner()[facei];
            const bool ownProtected = unrefineableCells.get(own);

            if
            (
                ownProtected
             && (neiLevel[facei-mesh().nInternalFaces()] > cellLevel[own])
            )
            {
                seedFace[facei] = true;
            }
        }

        syncTools::syncFaceList(mesh(), seedFace, orEqOp<bool>());


        // Extend unrefineableCells
        bool hasExtended = false;

        for (label facei = 0; facei < mesh().nInternalFaces(); facei++)
        {
            if (seedFace[facei])
            {
                const label own = mesh().faceOwner()[facei];
                if (unrefineableCells.get(own) == 0)
                {
                    unrefineableCells.set(own, 1);
                    hasExtended = true;
                }

                const label nei = mesh().faceNeighbour()[facei];
                if (unrefineableCells.get(nei) == 0)
                {
                    unrefineableCells.set(nei, 1);
                    hasExtended = true;
                }
            }
        }

        for
        (
            label facei = mesh().nInternalFaces();
            facei < mesh().nFaces();
            facei++
        )
        {
            if (seedFace[facei])
            {
                const label own = mesh().faceOwner()[facei];

                if (unrefineableCells.get(own) == 0)
                {
                    unrefineableCells.set(own, 1);
                    hasExtended = true;
                }
            }
        }

        if (!returnReduce(hasExtended, orOp<bool>()))
        {
            break;
        }
    }
}


void Foam::fvMeshTopoChangers::refiner::readDict()
{
    refineInterval_ = dict_.lookup<label>("refineInterval");

    if (refineInterval_ < 0)
    {
        FatalIOErrorInFunction(dict_)
            << "Illegal refineInterval " << refineInterval_ << nl
            << "The refineInterval setting in the dynamicMeshDict should"
            << " be >= 1." << nl
            << exit(FatalIOError);
    }

    maxCells_ = dict_.lookup<label>("maxCells");

    if (maxCells_ <= 0)
    {
        FatalIOErrorInFunction(dict_)
            << "Illegal maximum number of cells " << maxCells_ << nl
            << "The maxCells setting in the dynamicMeshDict should"
            << " be > 0." << nl
            << exit(FatalIOError);
    }

    nBufferLayers_ = dict_.lookup<label>("nBufferLayers");

    const List<Pair<word>> fluxVelocities = List<Pair<word>>
    (
        dict_.lookup("correctFluxes")
    );

    // Rework into hashtable.
    correctFluxes_.resize(fluxVelocities.size());
    forAll(fluxVelocities, i)
    {
        correctFluxes_.insert(fluxVelocities[i][0], fluxVelocities[i][1]);
    }

    dumpLevel_ = Switch(dict_.lookup("dumpLevel"));
}


Foam::autoPtr<Foam::polyTopoChangeMap>
Foam::fvMeshTopoChangers::refiner::refine
(
    const labelList& cellsToRefine
)
{
    // Mesh changing engine.
    polyTopoChange meshMod(mesh());

    // Play refinement commands into mesh changer.
    meshCutter_.setRefinement(cellsToRefine, meshMod);

    // Create mesh (with inflation), return map from old to new mesh.
    // autoPtr<polyTopoChangeMap> map = meshMod.changeMesh(*this, true);
    autoPtr<polyTopoChangeMap> map = meshMod.changeMesh(mesh(), false);

    Info<< "Refined from "
        << returnReduce(map().nOldCells(), sumOp<label>())
        << " to " << mesh().globalData().nTotalCells() << " cells." << endl;

    if (debug)
    {
        // Check map.
        for (label facei = 0; facei < mesh().nInternalFaces(); facei++)
        {
            const label oldFacei = map().faceMap()[facei];

            if (oldFacei >= mesh().nInternalFaces())
            {
                FatalErrorInFunction
                    << "New internal face:" << facei
                    << " fc:" << mesh().faceCentres()[facei]
                    << " originates from boundary oldFace:" << oldFacei
                    << abort(FatalError);
            }
        }
    }

    // Update fields
    mesh().topoChange(map);

    {
        // Correct the flux for modified/added faces. All the faces which only
        // have been renumbered will already have been handled by the mapping.
        const labelList& faceMap = map().faceMap();
        const labelList& reverseFaceMap = map().reverseFaceMap();

        // Storage for any master faces. These will be the original faces
        // on the coarse cell that get split into four (or rather the
        // master face gets modified and three faces get added from the master)
        labelHashSet masterFaces(4*cellsToRefine.size());

        forAll(faceMap, facei)
        {
            const label oldFacei = faceMap[facei];

            if (oldFacei >= 0)
            {
                const label masterFacei = reverseFaceMap[oldFacei];

                if (masterFacei < 0)
                {
                    FatalErrorInFunction
                        << "Problem: should not have removed faces"
                        << " when refining."
                        << nl << "face:" << facei << abort(FatalError);
                }
                else if (masterFacei != facei)
                {
                    masterFaces.insert(masterFacei);
                }
            }
        }
        if (debug)
        {
            Pout<< "Found " << masterFaces.size() << " split faces " << endl;
        }

        refineFluxes(masterFaces, map());
        refineUfs(masterFaces, map());
    }

    // Update numbering of protectedCells_
    if (protectedCells_.size())
    {
        PackedBoolList newProtectedCell(mesh().nCells());

        forAll(newProtectedCell, celli)
        {
            const label oldCelli = map().cellMap()[celli];
            newProtectedCell.set(celli, protectedCells_.get(oldCelli));
        }
        protectedCells_.transfer(newProtectedCell);
    }

    // Debug: Check refinement levels (across faces only)
    meshCutter_.checkRefinementLevels(-1, labelList(0));

    return map;
}


Foam::autoPtr<Foam::polyTopoChangeMap>
Foam::fvMeshTopoChangers::refiner::unrefine
(
    const labelList& splitPoints
)
{
    polyTopoChange meshMod(mesh());

    // Play refinement commands into mesh changer.
    meshCutter_.setUnrefinement(splitPoints, meshMod);


    // Save information on faces that will be combined
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Find the faceMidPoints on cells to be combined.
    // for each face resulting of split of face into four store the
    // midpoint
    Map<label> faceToSplitPoint(3*splitPoints.size());

    {
        forAll(splitPoints, i)
        {
            const label pointi = splitPoints[i];
            const labelList& pEdges = mesh().pointEdges()[pointi];

            forAll(pEdges, j)
            {
                const label otherPointi =
                    mesh().edges()[pEdges[j]].otherVertex(pointi);

                const labelList& pFaces = mesh().pointFaces()[otherPointi];

                forAll(pFaces, pFacei)
                {
                    faceToSplitPoint.insert(pFaces[pFacei], otherPointi);
                }
            }
        }
    }


    // Change mesh and generate map.
    // autoPtr<polyTopoChangeMap> map = meshMod.changeMesh(mesh(), true);
    autoPtr<polyTopoChangeMap> map = meshMod.changeMesh(mesh(), false);

    Info<< "Unrefined from "
        << returnReduce(map().nOldCells(), sumOp<label>())
        << " to " << mesh().globalData().nTotalCells() << " cells."
        << endl;

    // Update fields
    mesh().topoChange(map);

    // Correct the fluxes for modified faces
    unrefineFluxes(faceToSplitPoint, map());

    // Correct the face velocities for modified faces
    unrefineUfs(faceToSplitPoint, map());

    // Update numbering of protectedCells_
    if (protectedCells_.size())
    {
        PackedBoolList newProtectedCell(mesh().nCells());

        forAll(newProtectedCell, celli)
        {
            const label oldCelli = map().cellMap()[celli];
            if (oldCelli >= 0)
            {
                newProtectedCell.set(celli, protectedCells_.get(oldCelli));
            }
        }
        protectedCells_.transfer(newProtectedCell);
    }

    // Debug: Check refinement levels (across faces only)
    meshCutter_.checkRefinementLevels(-1, labelList(0));

    return map;
}


Foam::word Foam::fvMeshTopoChangers::refiner::Uname
(
    const surfaceVectorField& Uf
) const
{
    const word UfName(Uf.member());

    return
        IOobject::groupName
        (
            UfName.back() == 'f'
          ? word(UfName(UfName.size() - 1))
          : UfName.compare(UfName.size() - 3, 3, "f_0") == 0
            ? word(UfName(UfName.size() - 3) + "_0")
            : word::null,
            Uf.group()
        );
}


void Foam::fvMeshTopoChangers::refiner::refineFluxes
(
    const labelHashSet& masterFaces,
    const polyTopoChangeMap& map
)
{
    // Correct the flux for modified/added faces. All the faces which only
    // have been renumbered will already have been handled by the mapping.
    const labelList& faceMap = map.faceMap();
    const labelList& reverseFaceMap = map.reverseFaceMap();

    HashTable<surfaceScalarField*> fluxes
    (
        mesh().lookupClass<surfaceScalarField>()
    );
    forAllIter(HashTable<surfaceScalarField*>, fluxes, iter)
    {
        if (!correctFluxes_.found(iter.key()))
        {
            WarningInFunction
                << "Cannot find surfaceScalarField " << iter.key()
                << " in user-provided flux mapping table "
                << correctFluxes_ << endl
                << "    The flux mapping table is used to recreate the"
                << " flux on newly created faces." << endl
                << "    Either add the entry if it is a flux or use ("
                << iter.key() << " none) to suppress this warning."
                << endl;
            continue;
        }

        const word& Uname = correctFluxes_[iter.key()];

        if (Uname == "none")
        {
            continue;
        }

        if (Uname == "NaN")
        {
            Pout<< "Setting surfaceScalarField " << iter.key()
                << " to NaN" << endl;

            surfaceScalarField& phi = *iter();

            sigFpe::fillNan(phi.primitiveFieldRef());

            continue;
        }

        if (debug)
        {
            Pout<< "Mapping flux " << iter.key()
                << " using interpolated flux " << Uname
                << endl;
        }

        surfaceScalarField& phi = *iter();
        const surfaceScalarField phiU
        (
            fvc::interpolate(mesh().lookupObject<volVectorField>(Uname))
          & mesh().Sf()
        );

        // Recalculate new internal faces.
        for (label facei = 0; facei < mesh().nInternalFaces(); facei++)
        {
            const label oldFacei = faceMap[facei];

            if (oldFacei == -1)
            {
                // Inflated/appended
                phi[facei] = phiU[facei];
            }
            else if (reverseFaceMap[oldFacei] != facei)
            {
                // face-from-masterface
                phi[facei] = phiU[facei];
            }
        }

        // Recalculate new boundary faces.
        surfaceScalarField::Boundary& phiBf = phi.boundaryFieldRef();

        forAll(phiBf, patchi)
        {
            fvsPatchScalarField& patchPhi = phiBf[patchi];

            const fvsPatchScalarField& patchPhiU =
                phiU.boundaryField()[patchi];

            label facei = patchPhi.patch().start();

            forAll(patchPhi, i)
            {
                const label oldFacei = faceMap[facei];

                if (oldFacei == -1)
                {
                    // Inflated/appended
                    patchPhi[i] = patchPhiU[i];
                }
                else if (reverseFaceMap[oldFacei] != facei)
                {
                    // face-from-masterface
                    patchPhi[i] = patchPhiU[i];
                }

                facei++;
            }
        }

        // Update master faces
        forAllConstIter(labelHashSet, masterFaces, iter)
        {
            const label facei = iter.key();

            if (mesh().isInternalFace(facei))
            {
                phi[facei] = phiU[facei];
            }
            else
            {
                const label patchi =
                    mesh().boundaryMesh().whichPatch(facei);

                const label i =
                    facei - mesh().boundaryMesh()[patchi].start();

                phiBf[patchi][i] = phiU.boundaryField()[patchi][i];
            }
        }
    }
}


void Foam::fvMeshTopoChangers::refiner::refineUfs
(
    const labelHashSet& masterFaces,
    const polyTopoChangeMap& map
)
{
    const labelList& faceMap = map.faceMap();
    const labelList& reverseFaceMap = map.reverseFaceMap();

    // Interpolate U to Uf for added faces
    HashTable<surfaceVectorField*> Ufs
    (
        mesh().lookupClass<surfaceVectorField>()
    );
    forAllIter(HashTable<surfaceVectorField*>, Ufs, iter)
    {
        surfaceVectorField& Uf = *iter();

        const word Uname(this->Uname(Uf));

        if (Uname != word::null)
        {
            const surfaceVectorField UfU
            (
                fvc::interpolate
                (
                    mesh().lookupObject<volVectorField>(Uname)
                )
            );

            // Recalculate new internal faces.
            for (label facei = 0; facei < mesh().nInternalFaces(); facei++)
            {
                label oldFacei = faceMap[facei];

                if (oldFacei == -1)
                {
                    // Inflated/appended
                    Uf[facei] = UfU[facei];
                }
                else if (reverseFaceMap[oldFacei] != facei)
                {
                    // face-from-masterface
                    Uf[facei] = UfU[facei];
                }
            }

            // Recalculate new boundary faces.
            surfaceVectorField::Boundary& UfBf = Uf.boundaryFieldRef();
            forAll(UfBf, patchi)
            {
                fvsPatchVectorField& patchUf = UfBf[patchi];
                const fvsPatchVectorField& patchUfU =
                    UfU.boundaryField()[patchi];

                label facei = patchUf.patch().start();

                forAll(patchUf, i)
                {
                    label oldFacei = faceMap[facei];

                    if (oldFacei == -1)
                    {
                        // Inflated/appended
                        patchUf[i] = patchUfU[i];
                    }
                    else if (reverseFaceMap[oldFacei] != facei)
                    {
                        // face-from-masterface
                        patchUf[i] = patchUfU[i];
                    }

                    facei++;
                }
            }

            // Update master faces
            forAllConstIter(labelHashSet, masterFaces, iter)
            {
                label facei = iter.key();

                if (mesh().isInternalFace(facei))
                {
                    Uf[facei] = UfU[facei];
                }
                else
                {
                    const label patchi =
                        mesh().boundaryMesh().whichPatch(facei);
                    const label i =
                        facei - mesh().boundaryMesh()[patchi].start();

                    const fvsPatchVectorField& patchUfU =
                        UfU.boundaryField()[patchi];

                    fvsPatchVectorField& patchUf = UfBf[patchi];

                    patchUf[i] = patchUfU[i];
                }
            }
        }
    }
}


void Foam::fvMeshTopoChangers::refiner::unrefineFluxes
(
    const Map<label>& faceToSplitPoint,
    const polyTopoChangeMap& map
)
{
    const labelList& reversePointMap = map.reversePointMap();
    const labelList& reverseFaceMap = map.reverseFaceMap();

    HashTable<surfaceScalarField*> fluxes
    (
        mesh().lookupClass<surfaceScalarField>()
    );
    forAllIter(HashTable<surfaceScalarField*>, fluxes, iter)
    {
        if (!correctFluxes_.found(iter.key()))
        {
            WarningInFunction
                << "Cannot find surfaceScalarField " << iter.key()
                << " in user-provided flux mapping table "
                << correctFluxes_ << endl
                << "    The flux mapping table is used to recreate the"
                << " flux on newly created faces." << endl
                << "    Either add the entry if it is a flux or use ("
                << iter.key() << " none) to suppress this warning."
                << endl;
            continue;
        }

        const word& Uname = correctFluxes_[iter.key()];

        if (Uname == "none")
        {
            continue;
        }

        if (debug)
        {
            Info<< "Mapping flux " << iter.key()
                << " using interpolated flux " << Uname
                << endl;
        }

        surfaceScalarField& phi = *iter();
        surfaceScalarField::Boundary& phiBf = phi.boundaryFieldRef();

        const surfaceScalarField phiU
        (
            fvc::interpolate(mesh().lookupObject<volVectorField>(Uname))
          & mesh().Sf()
        );

        forAllConstIter(Map<label>, faceToSplitPoint, iter)
        {
            const label oldFacei = iter.key();
            const label oldPointi = iter();

            if (reversePointMap[oldPointi] < 0)
            {
                // midpoint was removed. See if face still exists.
                const label facei = reverseFaceMap[oldFacei];

                if (facei >= 0)
                {
                    if (mesh().isInternalFace(facei))
                    {
                        phi[facei] = phiU[facei];
                    }
                    else
                    {
                        const label patchi =
                            mesh().boundaryMesh().whichPatch(facei);

                        const label i =
                            facei - mesh().boundaryMesh()[patchi].start();

                        phiBf[patchi][i] = phiU.boundaryField()[patchi][i];
                    }
                }
            }
        }
    }
}


void Foam::fvMeshTopoChangers::refiner::unrefineUfs
(
    const Map<label>& faceToSplitPoint,
    const polyTopoChangeMap& map
)
{
    const labelList& reversePointMap = map.reversePointMap();
    const labelList& reverseFaceMap = map.reverseFaceMap();

    // Interpolate U to Uf for added faces
    HashTable<surfaceVectorField*> Ufs
    (
        mesh().lookupClass<surfaceVectorField>()
    );
    forAllIter(HashTable<surfaceVectorField*>, Ufs, iter)
    {
        surfaceVectorField& Uf = *iter();

        const word Uname(this->Uname(Uf));

        if (Uname != word::null)
        {
            surfaceVectorField::Boundary& UfBf = Uf.boundaryFieldRef();

            const surfaceVectorField UfU
            (
                fvc::interpolate(mesh().lookupObject<volVectorField>(Uname))
            );

            forAllConstIter(Map<label>, faceToSplitPoint, iter)
            {
                const label oldFacei = iter.key();
                const label oldPointi = iter();

                if (reversePointMap[oldPointi] < 0)
                {
                    // midpoint was removed. See if face still exists.
                    const label facei = reverseFaceMap[oldFacei];

                    if (facei >= 0)
                    {
                        if (mesh().isInternalFace(facei))
                        {
                            Uf[facei] = UfU[facei];
                        }
                        else
                        {
                            const label patchi =
                                mesh().boundaryMesh().whichPatch(facei);
                            const label i =
                                facei - mesh().boundaryMesh()[patchi].start();

                            UfBf[patchi][i] = UfU.boundaryField()[patchi][i];
                        }
                    }
                }
            }
        }
    }
}


const Foam::cellZone& Foam::fvMeshTopoChangers::refiner::findCellZone
(
    const word& cellZoneName
) const
{
    const label cellZoneID = mesh().cellZones().findZoneID(cellZoneName);

    bool cellZoneFound = (cellZoneID != -1);
    reduce(cellZoneFound, orOp<bool>());

    if (!cellZoneFound)
    {
        FatalErrorInFunction
            << "cannot find cellZone " << cellZoneName
            << exit(FatalError);
    }

    return mesh().cellZones()[cellZoneID];
}


Foam::scalarField
Foam::fvMeshTopoChangers::refiner::cellToPoint(const scalarField& vFld) const
{
    scalarField pFld(mesh().nPoints());

    forAll(mesh().pointCells(), pointi)
    {
        const labelList& pCells = mesh().pointCells()[pointi];

        scalar sum = 0.0;
        forAll(pCells, i)
        {
            sum += vFld[pCells[i]];
        }
        pFld[pointi] = sum/pCells.size();
    }

    return pFld;
}


Foam::scalarField Foam::fvMeshTopoChangers::refiner::error
(
    const scalarField& fld,
    const scalar minLevel,
    const scalar maxLevel
) const
{
    scalarField c(fld.size(), -1);

    forAll(c, celli)
    {
        scalar err = min(fld[celli] - minLevel, maxLevel - fld[celli]);

        if (err >= 0)
        {
            c[celli] = err;
        }
    }

    return c;
}


Foam::scalarField Foam::fvMeshTopoChangers::refiner::error
(
    const scalarField& fld,
    const labelList& cells,
    const scalar minLevel,
    const scalar maxLevel
) const
{
    scalarField c(fld.size(), -1);

    forAll(cells, i)
    {
        const label celli = cells[i];

        scalar err = min(fld[celli] - minLevel, maxLevel - fld[celli]);

        if (err >= 0)
        {
            c[celli] = err;
        }
    }

    return c;
}


void Foam::fvMeshTopoChangers::refiner::selectRefineCandidates
(
    PackedBoolList& candidateCells,
    const scalar lowerRefineLevel,
    const scalar upperRefineLevel,
    const scalar maxRefinement,
    const scalarField& vFld
) const
{
    // Get error per cell. Is -1 (not to be refined) to >0 (to be refined,
    // higher more desirable to be refined).
    const scalarField cellError
    (
        error(vFld, lowerRefineLevel, upperRefineLevel)
    );

    // Mark cells that are candidates for refinement.
    forAll(cellError, celli)
    {
        if (cellError[celli] > 0)
        {
            candidateCells.set(celli, 1);
        }
    }
}


void Foam::fvMeshTopoChangers::refiner::selectRefineCandidates
(
    PackedBoolList& candidateCells,
    const scalar lowerRefineLevel,
    const scalar upperRefineLevel,
    const scalar maxRefinement,
    const scalarField& vFld,
    const labelList& cells
) const
{
    // Get error per cell. Is -1 (not to be refined) to >0 (to be refined,
    // higher more desirable to be refined).
    const scalarField cellError
    (
        error(vFld, cells, lowerRefineLevel, upperRefineLevel)
    );

    // Mark cells that are candidates for refinement.
    forAll(cellError, celli)
    {
        if (cellError[celli] > 0)
        {
            candidateCells.set(celli, 1);
        }
    }
}


Foam::scalar Foam::fvMeshTopoChangers::refiner::selectRefineCandidates
(
    PackedBoolList& candidateCells,
    const dictionary& refineDict
) const
{
    const word fieldName(refineDict.lookup("field"));

    const volScalarField& vFld = mesh().lookupObject<volScalarField>(fieldName);

    const scalar lowerRefineLevel =
        refineDict.lookup<scalar>("lowerRefineLevel");
    const scalar upperRefineLevel =
        refineDict.lookup<scalar>("upperRefineLevel");

    const label maxRefinement = refineDict.lookup<label>("maxRefinement");

    if (maxRefinement <= 0)
    {
        FatalErrorInFunction
            << "Illegal maximum refinement level " << maxRefinement << nl
            << "The maxCells setting in the dynamicMeshDict should"
            << " be > 0." << nl
            << exit(FatalError);
    }

    if (refineDict.found("cellZone"))
    {
        // Determine candidates for refinement (looking at field only)
        selectRefineCandidates
        (
            candidateCells,
            lowerRefineLevel,
            upperRefineLevel,
            maxRefinement,
            vFld,
            findCellZone(refineDict.lookup("cellZone"))
        );
    }
    else
    {
        // Determine candidates for refinement (looking at field only)
        selectRefineCandidates
        (
            candidateCells,
            lowerRefineLevel,
            upperRefineLevel,
            maxRefinement,
            vFld
        );
    }

    return maxRefinement;
}


Foam::labelList Foam::fvMeshTopoChangers::refiner::selectRefineCells
(
    const label maxCells,
    const label maxRefinement,
    const PackedBoolList& candidateCells
) const
{
    // Every refined cell causes 7 extra cells
    const label nTotToRefine = (maxCells - mesh().globalData().nTotalCells())/7;

    const labelList& cellLevel = meshCutter_.cellLevel();

    // Mark cells that cannot be refined since they would trigger refinement
    // of protected cells (since 2:1 cascade)
    PackedBoolList unrefineableCells;
    calculateProtectedCells(unrefineableCells);

    // Count current selection
    const label nLocalCandidates = count(candidateCells, 1);
    const label nCandidates = returnReduce(nLocalCandidates, sumOp<label>());

    // Collect all cells
    DynamicList<label> candidates(nLocalCandidates);

    if (nCandidates < nTotToRefine)
    {
        forAll(candidateCells, celli)
        {
            if
            (
                candidateCells.get(celli)
             && (
                    unrefineableCells.empty()
                 || !unrefineableCells.get(celli)
                )
            )
            {
                candidates.append(celli);
            }
        }
    }
    else
    {
        // Sort by error? For now just truncate.
        for (label level = 0; level < maxRefinement; level++)
        {
            forAll(candidateCells, celli)
            {
                if
                (
                    cellLevel[celli] == level
                 && candidateCells.get(celli)
                 && (
                        unrefineableCells.empty()
                     || !unrefineableCells.get(celli)
                    )
                )
                {
                    candidates.append(celli);
                }
            }

            if (returnReduce(candidates.size(), sumOp<label>()) > nTotToRefine)
            {
                break;
            }
        }
    }

    // Guarantee 2:1 refinement after refinement
    labelList consistentSet
    (
        meshCutter_.consistentRefinement
        (
            candidates.shrink(),
            true               // Add to set to guarantee 2:1
        )
    );

    Info<< "Selected " << returnReduce(consistentSet.size(), sumOp<label>())
        << " cells for refinement out of " << mesh().globalData().nTotalCells()
        << "." << endl;

    return consistentSet;
}


Foam::labelList Foam::fvMeshTopoChangers::refiner::selectUnrefinePoints
(
    const PackedBoolList& markedCell
) const
{
    // All points that can be unrefined
    const labelList splitPoints(meshCutter_.getSplitPoints());

    DynamicList<label> newSplitPoints(splitPoints.size());

    forAll(splitPoints, i)
    {
        const label pointi = splitPoints[i];

        // Check that all cells are not marked
        const labelList& pCells = mesh().pointCells()[pointi];

        bool hasMarked = false;

        forAll(pCells, pCelli)
        {
            if (markedCell.get(pCells[pCelli]))
            {
                hasMarked = true;
                break;
            }
        }

        if (!hasMarked)
        {
            newSplitPoints.append(pointi);
        }
    }


    newSplitPoints.shrink();

    // Guarantee 2:1 refinement after unrefinement
    labelList consistentSet
    (
        meshCutter_.consistentUnrefinement
        (
            newSplitPoints,
            false
        )
    );

    Info<< "Selected " << returnReduce(consistentSet.size(), sumOp<label>())
        << " split points out of a possible "
        << returnReduce(splitPoints.size(), sumOp<label>())
        << "." << endl;

    return consistentSet;
}


void Foam::fvMeshTopoChangers::refiner::extendMarkedCells
(
    PackedBoolList& markedCell
) const
{
    // Mark faces using any marked cell
    boolList markedFace(mesh().nFaces(), false);

    forAll(markedCell, celli)
    {
        if (markedCell.get(celli))
        {
            const cell& cFaces = mesh().cells()[celli];

            forAll(cFaces, i)
            {
                markedFace[cFaces[i]] = true;
            }
        }
    }

    syncTools::syncFaceList(mesh(), markedFace, orEqOp<bool>());

    // Update cells using any markedFace
    for (label facei = 0; facei < mesh().nInternalFaces(); facei++)
    {
        if (markedFace[facei])
        {
            markedCell.set(mesh().faceOwner()[facei], 1);
            markedCell.set(mesh().faceNeighbour()[facei], 1);
        }
    }

    for
    (
        label facei = mesh().nInternalFaces();
        facei < mesh().nFaces();
        facei++
    )
    {
        if (markedFace[facei])
        {
            markedCell.set(mesh().faceOwner()[facei], 1);
        }
    }
}


void Foam::fvMeshTopoChangers::refiner::checkEightAnchorPoints
(
    PackedBoolList& protectedCell,
    label& nProtected
) const
{
    const labelList& cellLevel = meshCutter_.cellLevel();
    const labelList& pointLevel = meshCutter_.pointLevel();

    labelList nAnchorPoints(mesh().nCells(), 0);

    forAll(pointLevel, pointi)
    {
        const labelList& pCells = mesh().pointCells(pointi);

        forAll(pCells, pCelli)
        {
            const label celli = pCells[pCelli];

            if (pointLevel[pointi] <= cellLevel[celli])
            {
                // Check if cell has already 8 anchor points -> protect cell
                if (nAnchorPoints[celli] == 8)
                {
                    if (protectedCell.set(celli, true))
                    {
                        nProtected++;
                    }
                }

                if (!protectedCell[celli])
                {
                    nAnchorPoints[celli]++;
                }
            }
        }
    }

    forAll(protectedCell, celli)
    {
        if (!protectedCell[celli] && nAnchorPoints[celli] != 8)
        {
            protectedCell.set(celli, true);
            nProtected++;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshTopoChangers::refiner::refiner(fvMesh& mesh, const dictionary& dict)
:
    fvMeshTopoChanger(mesh),
    dict_(dict),
    meshCutter_(mesh),
    dumpLevel_(false),
    nRefinementIterations_(0),
    protectedCells_(mesh.nCells(), 0),
    changedSinceWrite_(false),
    timeIndex_(-1)
{
    // Read static part of dictionary
    readDict();

    const labelList& cellLevel = meshCutter_.cellLevel();
    const labelList& pointLevel = meshCutter_.pointLevel();

    // Set cells that should not be refined.
    // This is currently any cell which does not have 8 anchor points or
    // uses any face which does not have 4 anchor points.
    // Note: do not use cellPoint addressing

    // Count number of points <= cellLevel
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList nAnchors(mesh.nCells(), 0);

    label nProtected = 0;

    forAll(mesh.pointCells(), pointi)
    {
        const labelList& pCells = mesh.pointCells()[pointi];

        forAll(pCells, i)
        {
            const label celli = pCells[i];

            if (!protectedCells_.get(celli))
            {
                if (pointLevel[pointi] <= cellLevel[celli])
                {
                    nAnchors[celli]++;

                    if (nAnchors[celli] > 8)
                    {
                        protectedCells_.set(celli, 1);
                        nProtected++;
                    }
                }
            }
        }
    }


    // Count number of points <= faceLevel
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Bit tricky since proc face might be one more refined than the owner since
    // the coupled one is refined.

    {
        labelList neiLevel(mesh.nFaces());

        for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
        {
            neiLevel[facei] = cellLevel[mesh.faceNeighbour()[facei]];
        }

        for
        (
            label facei = mesh.nInternalFaces();
            facei < mesh.nFaces();
            facei++
        )
        {
            neiLevel[facei] = cellLevel[mesh.faceOwner()[facei]];
        }
        syncTools::swapFaceList(mesh, neiLevel);


        boolList protectedFace(mesh.nFaces(), false);

        forAll(mesh.faceOwner(), facei)
        {
            const label faceLevel = max
            (
                cellLevel[mesh.faceOwner()[facei]],
                neiLevel[facei]
            );

            const face& f = mesh.faces()[facei];

            label nAnchors = 0;

            forAll(f, fp)
            {
                if (pointLevel[f[fp]] <= faceLevel)
                {
                    nAnchors++;

                    if (nAnchors > 4)
                    {
                        protectedFace[facei] = true;
                        break;
                    }
                }
            }
        }

        syncTools::syncFaceList(mesh, protectedFace, orEqOp<bool>());

        for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
        {
            if (protectedFace[facei])
            {
                protectedCells_.set(mesh.faceOwner()[facei], 1);
                nProtected++;
                protectedCells_.set(mesh.faceNeighbour()[facei], 1);
                nProtected++;
            }
        }

        for
        (
            label facei = mesh.nInternalFaces();
            facei < mesh.nFaces();
            facei++
        )
        {
            if (protectedFace[facei])
            {
                protectedCells_.set(mesh.faceOwner()[facei], 1);
                nProtected++;
            }
        }

        // Also protect any cells that are less than hex
        forAll(mesh.cells(), celli)
        {
            const cell& cFaces = mesh.cells()[celli];

            if (cFaces.size() < 6)
            {
                if (protectedCells_.set(celli, 1))
                {
                    nProtected++;
                }
            }
            else
            {
                forAll(cFaces, cFacei)
                {
                    if (mesh.faces()[cFaces[cFacei]].size() < 4)
                    {
                        if (protectedCells_.set(celli, 1))
                        {
                            nProtected++;
                        }
                        break;
                    }
                }
            }
        }

        // Check cells for 8 corner points
        checkEightAnchorPoints(protectedCells_, nProtected);
    }

    if (returnReduce(nProtected, sumOp<label>()) == 0)
    {
        protectedCells_.clear();
    }
    else
    {
        cellSet protectedCells(mesh, "protectedCells", nProtected);
        forAll(protectedCells_, celli)
        {
            if (protectedCells_[celli])
            {
                protectedCells.insert(celli);
            }
        }

        Info<< "Detected " << returnReduce(nProtected, sumOp<label>())
            << " cells that are protected from refinement."
            << " Writing these to cellSet "
            << protectedCells.name()
            << "." << endl;

        protectedCells.write();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshTopoChangers::refiner::~refiner()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fvMeshTopoChangers::refiner::update()
{
    // Only refine on the first call in a time-step
    if (timeIndex_ != mesh().time().timeIndex())
    {
        timeIndex_ = mesh().time().timeIndex();
    }
    else
    {
        return false;
    }

    bool hasChanged = false;

    if (refineInterval_ == 0)
    {
        return hasChanged;
    }

    // Note: cannot refine at time 0 since no V0 present since mesh not
    //       moved yet.

    if
    (
        mesh().time().timeIndex() > 0
     && mesh().time().timeIndex() % refineInterval_ == 0
    )
    {
        // Cells marked for refinement or otherwise protected from unrefinement.
        PackedBoolList refineCells(mesh().nCells());

        label maxRefinement = 0;

        if (dict_.isDict("refinementRegions"))
        {
            const dictionary& refinementRegions
            (
                dict_.subDict("refinementRegions")
            );

            forAllConstIter(dictionary, refinementRegions, iter)
            {
                maxRefinement = max
                (
                    selectRefineCandidates
                    (
                        refineCells,
                        refinementRegions.subDict(iter().keyword())
                    ),
                    maxRefinement
                );
            }
        }
        else
        {
            maxRefinement = selectRefineCandidates(refineCells, dict_);
        }

        // Extend with a buffer layer to prevent neighbouring points
        // being unrefined.
        for (label i = 0; i < nBufferLayers_; i++)
        {
            extendMarkedCells(refineCells);
        }

        PackedBoolList refinableCells(refineCells);

        {
            const labelList& cellLevel = meshCutter_.cellLevel();

            // Mark cells that are candidates for refinement.
            forAll(cellLevel, celli)
            {
                if (cellLevel[celli] >= maxRefinement)
                {
                    refinableCells.unset(celli);
                }
            }
        }

        if (mesh().globalData().nTotalCells() < maxCells_)
        {
            // Select subset of candidates. Take into account max allowable
            // cells, refinement level, protected cells.
            const labelList cellsToRefine
            (
                selectRefineCells
                (
                    maxCells_,
                    maxRefinement,
                    refinableCells
                )
            );

            const label nCellsToRefine = returnReduce
            (
                cellsToRefine.size(), sumOp<label>()
            );

            if (nCellsToRefine > 0)
            {
                // Refine/update mesh and map fields
                autoPtr<polyTopoChangeMap> map = refine(cellsToRefine);

                // Update refinableCells. Note that some of the marked ones have
                // not been refined due to constraints.
                {
                    const labelList& cellMap = map().cellMap();
                    const labelList& reverseCellMap = map().reverseCellMap();

                    PackedBoolList newRefineCell(cellMap.size());

                    forAll(cellMap, celli)
                    {
                        const label oldCelli = cellMap[celli];

                        if (oldCelli < 0)
                        {
                            newRefineCell.set(celli, 1);
                        }
                        else if (reverseCellMap[oldCelli] != celli)
                        {
                            newRefineCell.set(celli, 1);
                        }
                        else
                        {
                            newRefineCell.set
                            (
                                celli,
                                refinableCells.get(oldCelli)
                            );
                        }
                    }
                    refinableCells.transfer(newRefineCell);
                }

                hasChanged = true;
            }
        }

        {
            // Select unrefineable points that are not marked in refineCells
            const labelList pointsToUnrefine(selectUnrefinePoints(refineCells));

            const label nSplitPoints = returnReduce
            (
                pointsToUnrefine.size(),
                sumOp<label>()
            );

            if (nSplitPoints > 0)
            {
                // Refine/update mesh
                unrefine(pointsToUnrefine);

                hasChanged = true;
            }
        }


        if ((nRefinementIterations_ % 10) == 0)
        {
            // Compact refinement history occasionally (how often?).
            // Unrefinement causes holes in the refinementHistory.
            const_cast<refinementHistory&>(meshCutter().history()).compact();
        }
        nRefinementIterations_++;
    }

    if (hasChanged)
    {
        changedSinceWrite_ = true;
    }

    return hasChanged;
}


void Foam::fvMeshTopoChangers::refiner::topoChange(const polyTopoChangeMap& map)
{
    // Update numbering of cells/vertices.
    meshCutter_.topoChange(map);
}


void Foam::fvMeshTopoChangers::refiner::mapMesh(const polyMeshMap& map)
{
    // meshCutter_ will need to be re-constructed from the new mesh
    // and protectedCells_ updated.
    // The constructor should be refactored for the protectedCells_ update.
    NotImplemented;
}


void Foam::fvMeshTopoChangers::refiner::distribute
(
    const polyDistributionMap& map
)
{
    // Redistribute the mesh cutting engine
    meshCutter_.distribute(map);
}


bool Foam::fvMeshTopoChangers::refiner::write(const bool write) const
{
    if (changedSinceWrite_)
    {
        // Force refinement data to go to the current time directory.
        const_cast<hexRef8&>(meshCutter_).setInstance(mesh().time().name());

        bool writeOk = meshCutter_.write(write);

        if (dumpLevel_)
        {
            volScalarField scalarCellLevel
            (
                IOobject
                (
                    "cellLevel",
                    mesh().time().name(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE,
                    false
                ),
                mesh(),
                dimensionedScalar(dimless, 0)
            );

            const labelList& cellLevel = meshCutter_.cellLevel();

            forAll(cellLevel, celli)
            {
                scalarCellLevel[celli] = cellLevel[celli];
            }

            writeOk = writeOk && scalarCellLevel.write();
        }

        changedSinceWrite_ = false;

        return writeOk;
    }
    else
    {
        return true;
    }
}


// ************************************************************************* //
