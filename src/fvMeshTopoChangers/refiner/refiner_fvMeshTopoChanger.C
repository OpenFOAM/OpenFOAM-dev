/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "refiner_fvMeshTopoChanger.H"
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


void Foam::fvMeshTopoChangers::refiner::calcProtectedCells
(
    PackedBoolList& protectedCells
) const
{
    const label nFaces = mesh().nFaces();
    const label nInternalFaces = mesh().nInternalFaces();

    const labelList& cellLevel = meshCutter_.cellLevel();
    const labelList& pointLevel = meshCutter_.pointLevel();

    // Determine the number of parent anchors (i.e., pointLevel < cellLevel)
    // and anchors (pointLevel <= cellLevel) in each cell
    labelList cellNParentAnchors(mesh().nCells(), 0);
    labelList cellNAnchors(mesh().nCells(), 0);
    forAll(mesh().pointCells(), pointi)
    {
        const labelList& pCells = mesh().pointCells()[pointi];
        forAll(pCells, pci)
        {
            const label celli = pCells[pci];
            if (pointLevel[pointi] < cellLevel[celli])
            {
                cellNParentAnchors[celli] ++;
            }
            if (pointLevel[pointi] <= cellLevel[celli])
            {
                cellNAnchors[celli] ++;
            }
        }
    }

    // A cell is (probably) refine-able if it has 8 anchors and 1 parent
    // anchor, unless this is an original cell (i.e., cellLevel == 0), in which
    // case it must have 0 parent anchors
    forAll(mesh().cells(), celli)
    {
        protectedCells[celli] =
            (cellNAnchors[celli] != 8)
         || (cellLevel[celli] > 0 && cellNParentAnchors[celli] != 1)
         || (cellLevel[celli] == 0 && cellNParentAnchors[celli] != 0);
    }

    // Construct a list of face-neighbour levels. Requires communication.
    labelList neiLevel(nFaces);
    for (label facei = 0; facei < nInternalFaces; facei ++)
    {
        neiLevel[facei] = cellLevel[mesh().faceNeighbour()[facei]];
    }
    for (label facei = nInternalFaces; facei < nFaces; facei ++)
    {
        neiLevel[facei] = cellLevel[mesh().faceOwner()[facei]];
    }
    syncTools::swapFaceList(mesh(), neiLevel);

    // Count the number of parent anchors and anchors in each face
    labelList faceNParentAnchors(nFaces, 0);
    labelList faceNAnchors(nFaces, 0);
    forAll(mesh().faces(), facei)
    {
        const label fLevel =
            max(cellLevel[mesh()().faceOwner()[facei]], neiLevel[facei]);
        const face& f = mesh().faces()[facei];
        forAll(f, fpi)
        {
            if (pointLevel[f[fpi]] < fLevel)
            {
                faceNParentAnchors[facei] ++;
            }
            if (pointLevel[f[fpi]] <= fLevel)
            {
                faceNAnchors[facei] ++;
            }
        }
    }

    // A face is (probably) refine-able if it has 4 anchors and no more than 1
    // parent anchor (*), unless this is an original face (i.e., faceLevel ==
    // 0), in which case it must have 0 parent anchors.
    //
    // (*) Note this is slightly different from the cell criteria. Faces
    // constructed inside cells have no parent anchors.
    //
    boolList faceProtected(nFaces, false);
    forAll(mesh().faces(), facei)
    {
        const label fLevel =
            max(cellLevel[mesh().faceOwner()[facei]], neiLevel[facei]);

        faceProtected[facei] =
            (faceNAnchors[facei] != 4)
         || (fLevel > 0 && faceNParentAnchors[facei] > 1)
         || (fLevel == 0 && faceNParentAnchors[facei] != 0);
    }

    // Another criteria... Points are created when adjacent cells are refined.
    // This means then should always be created along an edge in predictable
    // patterns (e.g., with levels 1 4 3 4 2 4 3 4 1). No two adjacent points
    // of a face should therefore have the same level, unless the edge is
    // unrefined (e.g., 1 1). So, if we find two such adjacent points, then we
    // also need to consider the face protected.
    forAll(mesh().faces(), facei)
    {
        const face& f = mesh().faces()[facei];
        const label fLevel =
            max(cellLevel[mesh().faceOwner()[facei]], neiLevel[facei]);
        forAll(f, fpi)
        {
            if
            (
                pointLevel[f[fpi]] > fLevel
             && pointLevel[f[fpi]] == pointLevel[f[f.fcIndex(fpi)]]
            )
            {
                faceProtected[facei] = true;
            }
        }
    }

    syncTools::syncFaceList(mesh(), faceProtected, orEqOp<bool>());

    // Any cell which has a protected face is also protected
    for (label facei = 0; facei < nFaces; facei ++)
    {
        const label own = mesh().faceOwner()[facei];
        protectedCells[own] = protectedCells[own] || faceProtected[facei];
        if (facei < nInternalFaces)
        {
            const label nei = mesh().faceNeighbour()[facei];
            protectedCells[nei] = protectedCells[nei] || faceProtected[facei];
        }
    }
}


void Foam::fvMeshTopoChangers::refiner::calcAdditionallyProtectedCells
(
    const PackedBoolList& protectedCells,
    PackedBoolList& additionallyProtectedCells
) const
{
    if (protectedCells_.empty())
    {
        additionallyProtectedCells.clear();
        return;
    }

    const labelList& cellLevel = meshCutter_.cellLevel();

    additionallyProtectedCells = protectedCells_;

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
            const bool ownProtected = additionallyProtectedCells.get(own);
            const label nei = mesh().faceNeighbour()[facei];
            const bool neiProtected = additionallyProtectedCells.get(nei);

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
            const bool ownProtected = additionallyProtectedCells.get(own);

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


        // Extend additionallyProtectedCells
        bool hasExtended = false;

        for (label facei = 0; facei < mesh().nInternalFaces(); facei++)
        {
            if (seedFace[facei])
            {
                const label own = mesh().faceOwner()[facei];
                if (additionallyProtectedCells.get(own) == 0)
                {
                    additionallyProtectedCells.set(own, 1);
                    hasExtended = true;
                }

                const label nei = mesh().faceNeighbour()[facei];
                if (additionallyProtectedCells.get(nei) == 0)
                {
                    additionallyProtectedCells.set(nei, 1);
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

                if (additionallyProtectedCells.get(own) == 0)
                {
                    additionallyProtectedCells.set(own, 1);
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

    if (dict_.found("correctFluxes"))
    {
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
    }

    dumpLevel_ = Switch(dict_.lookup("dumpLevel"));
}


Foam::autoPtr<Foam::polyTopoChangeMap>
Foam::fvMeshTopoChangers::refiner::refine
(
    const labelList& cellsToRefine
)
{
    mesh().preChange();

    // Mesh changing engine.
    polyTopoChange meshMod(mesh());

    // Play refinement commands into mesh changer.
    meshCutter_.setRefinement(cellsToRefine, meshMod);

    // Create mesh
    // return map from old to new mesh.
    autoPtr<polyTopoChangeMap> map = meshMod.changeMesh(mesh());

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
    mesh().preChange();

    // Mesh changing engine.
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
    autoPtr<polyTopoChangeMap> map = meshMod.changeMesh(mesh());

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
    if (correctFluxes_.size())
    {
        UPtrList<surfaceScalarField> fluxes
        (
            mesh().fields<surfaceScalarField>()
        );

        forAll(fluxes, i)
        {
            surfaceScalarField& flux = fluxes[i];

            if (!correctFluxes_.found(flux.name()))
            {
                WarningInFunction
                    << "Cannot find surfaceScalarField " << flux.name()
                    << " in user-provided flux mapping table "
                    << correctFluxes_ << endl
                    << "    The flux mapping table is used to recreate the"
                    << " flux on newly created faces." << endl
                    << "    Either add the entry if it is a flux or use ("
                    << flux.name() << " none) to suppress this warning."
                    << endl;
            }
            else
            {
                const word& method = correctFluxes_[flux.name()];

                if (method == "none")
                {}
                else if (method == "NaN")
                {
                    Pout<< "Setting surfaceScalarField " << flux.name()
                        << " to NaN" << endl;

                    sigFpe::fillNan(flux.primitiveFieldRef());
                }
                else
                {
                    FatalErrorInFunction
                        << "Unknown refinement method " << method
                        << " for surfaceScalarField " << flux.name()
                        << " in user-provided flux mapping table "
                        << correctFluxes_
                        << exit(FatalError);
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
    if (correctFluxes_.size())
    {
        UPtrList<surfaceScalarField> fluxes
        (
            mesh().fields<surfaceScalarField>()
        );

        forAll(fluxes, i)
        {
            surfaceScalarField& flux = fluxes[i];

            if (!correctFluxes_.found(flux.name()))
            {
                WarningInFunction
                    << "Cannot find surfaceScalarField " << flux.name()
                    << " in user-provided flux mapping table "
                    << correctFluxes_ << endl
                    << "    The flux mapping table is used to recreate the"
                    << " flux on newly created faces." << endl
                    << "    Either add the entry if it is a flux or use ("
                    << flux.name() << " none) to suppress this warning."
                    << endl;
            }
            else
            {
                const word& method = correctFluxes_[flux.name()];

                if (method != "none")
                {
                    FatalErrorInFunction
                        << "Unknown unrefinement method " << method
                        << " for surfaceScalarField " << flux.name()
                        << " in user-provided flux mapping table "
                        << correctFluxes_
                        << exit(FatalError);
                }
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
    UPtrList<surfaceVectorField> Ufs(mesh().fields<surfaceVectorField>());

    forAll(Ufs, i)
    {
        surfaceVectorField& Uf = Ufs[i];

        const word Uname(this->Uname(Uf));

        if (Uname != word::null)
        {
            const surfaceVectorField UfU
            (
                fvc::interpolate(mesh().lookupObject<volVectorField>(Uname))
            );

            // Recalculate new internal faces.
            for (label facei = 0; facei < mesh().nInternalFaces(); facei++)
            {
                label oldFacei = faceMap[facei];

                if (oldFacei == -1)
                {
                    // Inserted
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
                        // Inserted/appended
                        patchUf[i] = patchUfU[i];
                    }
                    else if (reverseFaceMap[oldFacei] != facei)
                    {
                        // face-from-master-face
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
                        mesh().poly().boundary().whichPatch(facei);
                    const label i =
                        facei - mesh().poly().boundary()[patchi].start();

                    const fvsPatchVectorField& patchUfU =
                        UfU.boundaryField()[patchi];

                    fvsPatchVectorField& patchUf = UfBf[patchi];

                    patchUf[i] = patchUfU[i];
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
    UPtrList<surfaceVectorField> Ufs(mesh().fields<surfaceVectorField>());

    forAll(Ufs, i)
    {
        surfaceVectorField& Uf = Ufs[i];

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
                                mesh().poly().boundary().whichPatch(facei);
                            const label i =
                                facei
                              - mesh().poly().boundary()[patchi].start();

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
    const label cellZoneID = mesh().cellZones().findIndex(cellZoneName);

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

    // Calculate cells that are additionally protected, as their refinement
    // would trigger refinement of protectedCells due to the 2:1 limit
    PackedBoolList additionallyProtectedCells;
    calcAdditionallyProtectedCells(protectedCells_, additionallyProtectedCells);

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
                    additionallyProtectedCells.empty()
                 || !additionallyProtectedCells.get(celli)
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
                        additionallyProtectedCells.empty()
                     || !additionallyProtectedCells.get(celli)
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshTopoChangers::refiner::refiner(fvMesh& mesh, const dictionary& dict)
:
    fvMeshTopoChanger(mesh),
    dict_(dict),
    meshCutter_(hexRef8::New(mesh)),
    dumpLevel_(false),
    nRefinementIterations_(0),
    protectedCells_(mesh.nCells(), 0),
    timeIndex_(-1)
{
    // Read static part of dictionary
    readDict();

    // Calculate and report the protected cells
    calcProtectedCells(protectedCells_);

    // Report
    labelList protectedCellsSet = protectedCells_.used();
    const label nProtectedCells =
        returnReduce(protectedCellsSet.size(), sumOp<label>());
    if (nProtectedCells)
    {
        Info<< "Detected " << nProtectedCells << " cells that are protected "
            << "from refinement. Writing these to cell-set 'protectedCells'."
            << endl;

        cellSet(mesh, "protectedCells", protectedCellsSet).write();
    }
    else
    {
        protectedCells_.clear();
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
            const_cast<refinementHistory&>(meshCutter_.history()).compact();
        }
        nRefinementIterations_++;
    }

    return hasChanged;
}


void Foam::fvMeshTopoChangers::refiner::topoChange(const polyTopoChangeMap& map)
{}


void Foam::fvMeshTopoChangers::refiner::mapMesh(const polyMeshMap& map)
{}


void Foam::fvMeshTopoChangers::refiner::distribute
(
    const polyDistributionMap& map
)
{}


// ************************************************************************* //
