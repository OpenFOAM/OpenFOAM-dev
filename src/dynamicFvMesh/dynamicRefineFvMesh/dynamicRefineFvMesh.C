/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "dynamicRefineFvMesh.H"
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
    defineTypeNameAndDebug(dynamicRefineFvMesh, 0);
    addToRunTimeSelectionTable(dynamicFvMesh, dynamicRefineFvMesh, IOobject);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::label Foam::dynamicRefineFvMesh::count
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


void Foam::dynamicRefineFvMesh::calculateProtectedCells
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
    labelList neiLevel(nFaces()-nInternalFaces());

    for (label facei = nInternalFaces(); facei < nFaces(); facei++)
    {
        neiLevel[facei-nInternalFaces()] = cellLevel[faceOwner()[facei]];
    }
    syncTools::swapBoundaryFaceList(*this, neiLevel);


    while (true)
    {
        // Pick up faces on border of protected cells
        boolList seedFace(nFaces(), false);

        forAll(faceNeighbour(), facei)
        {
            label own = faceOwner()[facei];
            bool ownProtected = unrefineableCells.get(own);
            label nei = faceNeighbour()[facei];
            bool neiProtected = unrefineableCells.get(nei);

            if (ownProtected && (cellLevel[nei] > cellLevel[own]))
            {
                seedFace[facei] = true;
            }
            else if (neiProtected && (cellLevel[own] > cellLevel[nei]))
            {
                seedFace[facei] = true;
            }
        }
        for (label facei = nInternalFaces(); facei < nFaces(); facei++)
        {
            label own = faceOwner()[facei];
            bool ownProtected = unrefineableCells.get(own);
            if
            (
                ownProtected
             && (neiLevel[facei-nInternalFaces()] > cellLevel[own])
            )
            {
                seedFace[facei] = true;
            }
        }

        syncTools::syncFaceList(*this, seedFace, orEqOp<bool>());


        // Extend unrefineableCells
        bool hasExtended = false;

        for (label facei = 0; facei < nInternalFaces(); facei++)
        {
            if (seedFace[facei])
            {
                label own = faceOwner()[facei];
                if (unrefineableCells.get(own) == 0)
                {
                    unrefineableCells.set(own, 1);
                    hasExtended = true;
                }

                label nei = faceNeighbour()[facei];
                if (unrefineableCells.get(nei) == 0)
                {
                    unrefineableCells.set(nei, 1);
                    hasExtended = true;
                }
            }
        }
        for (label facei = nInternalFaces(); facei < nFaces(); facei++)
        {
            if (seedFace[facei])
            {
                label own = faceOwner()[facei];
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


void Foam::dynamicRefineFvMesh::readDict()
{
    const dictionary refineDict
    (
        dynamicMeshDict().optionalSubDict(typeName + "Coeffs")
    );

    List<Pair<word>> fluxVelocities = List<Pair<word>>
    (
        refineDict.lookup("correctFluxes")
    );
    // Rework into hashtable.
    correctFluxes_.resize(fluxVelocities.size());
    forAll(fluxVelocities, i)
    {
        correctFluxes_.insert(fluxVelocities[i][0], fluxVelocities[i][1]);
    }

    dumpLevel_ = Switch(refineDict.lookup("dumpLevel"));
}


// Refines cells, maps fields and recalculates (an approximate) flux
Foam::autoPtr<Foam::mapPolyMesh>
Foam::dynamicRefineFvMesh::refine
(
    const labelList& cellsToRefine
)
{
    // Mesh changing engine.
    polyTopoChange meshMod(*this);

    // Play refinement commands into mesh changer.
    meshCutter_.setRefinement(cellsToRefine, meshMod);

    // Create mesh (with inflation), return map from old to new mesh.
    // autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, true);
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, false);

    Info<< "Refined from "
        << returnReduce(map().nOldCells(), sumOp<label>())
        << " to " << globalData().nTotalCells() << " cells." << endl;

    if (debug)
    {
        // Check map.
        for (label facei = 0; facei < nInternalFaces(); facei++)
        {
            label oldFacei = map().faceMap()[facei];

            if (oldFacei >= nInternalFaces())
            {
                FatalErrorInFunction
                    << "New internal face:" << facei
                    << " fc:" << faceCentres()[facei]
                    << " originates from boundary oldFace:" << oldFacei
                    << abort(FatalError);
            }
        }
    }

    // Update fields
    updateMesh(map);

    // Correct the flux for modified/added faces. All the faces which only
    // have been renumbered will already have been handled by the mapping.
    {
        const labelList& faceMap = map().faceMap();
        const labelList& reverseFaceMap = map().reverseFaceMap();

        // Storage for any master faces. These will be the original faces
        // on the coarse cell that get split into four (or rather the
        // master face gets modified and three faces get added from the master)
        labelHashSet masterFaces(4*cellsToRefine.size());

        forAll(faceMap, facei)
        {
            label oldFacei = faceMap[facei];

            if (oldFacei >= 0)
            {
                label masterFacei = reverseFaceMap[oldFacei];

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

        HashTable<surfaceScalarField*> fluxes
        (
            lookupClass<surfaceScalarField>()
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

            const word& UName = correctFluxes_[iter.key()];

            if (UName == "none")
            {
                continue;
            }

            if (UName == "NaN")
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
                    << " using interpolated flux " << UName
                    << endl;
            }

            surfaceScalarField& phi = *iter();
            const surfaceScalarField phiU
            (
                fvc::interpolate
                (
                    lookupObject<volVectorField>(UName)
                )
              & Sf()
            );

            // Recalculate new internal faces.
            for (label facei = 0; facei < nInternalFaces(); facei++)
            {
                label oldFacei = faceMap[facei];

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
            surfaceScalarField::Boundary& phiBf =
                phi.boundaryFieldRef();
            forAll(phiBf, patchi)
            {
                fvsPatchScalarField& patchPhi = phiBf[patchi];
                const fvsPatchScalarField& patchPhiU =
                    phiU.boundaryField()[patchi];

                label facei = patchPhi.patch().start();

                forAll(patchPhi, i)
                {
                    label oldFacei = faceMap[facei];

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
                label facei = iter.key();

                if (isInternalFace(facei))
                {
                    phi[facei] = phiU[facei];
                }
                else
                {
                    label patchi = boundaryMesh().whichPatch(facei);
                    label i = facei - boundaryMesh()[patchi].start();

                    const fvsPatchScalarField& patchPhiU =
                        phiU.boundaryField()[patchi];

                    fvsPatchScalarField& patchPhi = phiBf[patchi];

                    patchPhi[i] = patchPhiU[i];
                }
            }
        }
    }


    // Update numbering of cells/vertices.
    meshCutter_.updateMesh(map);

    // Update numbering of protectedCells_
    if (protectedCells_.size())
    {
        PackedBoolList newProtectedCell(nCells());

        forAll(newProtectedCell, celli)
        {
            label oldCelli = map().cellMap()[celli];
            newProtectedCell.set(celli, protectedCells_.get(oldCelli));
        }
        protectedCells_.transfer(newProtectedCell);
    }

    // Debug: Check refinement levels (across faces only)
    meshCutter_.checkRefinementLevels(-1, labelList(0));

    return map;
}


Foam::autoPtr<Foam::mapPolyMesh>
Foam::dynamicRefineFvMesh::unrefine
(
    const labelList& splitPoints
)
{
    polyTopoChange meshMod(*this);

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
            label pointi = splitPoints[i];

            const labelList& pEdges = pointEdges()[pointi];

            forAll(pEdges, j)
            {
                label otherPointi = edges()[pEdges[j]].otherVertex(pointi);

                const labelList& pFaces = pointFaces()[otherPointi];

                forAll(pFaces, pFacei)
                {
                    faceToSplitPoint.insert(pFaces[pFacei], otherPointi);
                }
            }
        }
    }


    // Change mesh and generate map.
    // autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, true);
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, false);

    Info<< "Unrefined from "
        << returnReduce(map().nOldCells(), sumOp<label>())
        << " to " << globalData().nTotalCells() << " cells."
        << endl;

    // Update fields
    updateMesh(map);

    // Correct the flux for modified faces.
    {
        const labelList& reversePointMap = map().reversePointMap();
        const labelList& reverseFaceMap = map().reverseFaceMap();

        HashTable<surfaceScalarField*> fluxes
        (
            lookupClass<surfaceScalarField>()
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

            const word& UName = correctFluxes_[iter.key()];

            if (UName == "none")
            {
                continue;
            }

            if (debug)
            {
                Info<< "Mapping flux " << iter.key()
                    << " using interpolated flux " << UName
                    << endl;
            }

            surfaceScalarField& phi = *iter();
            surfaceScalarField::Boundary& phiBf =
                phi.boundaryFieldRef();

            const surfaceScalarField phiU
            (
                fvc::interpolate
                (
                    lookupObject<volVectorField>(UName)
                )
              & Sf()
            );


            forAllConstIter(Map<label>, faceToSplitPoint, iter)
            {
                label oldFacei = iter.key();
                label oldPointi = iter();

                if (reversePointMap[oldPointi] < 0)
                {
                    // midpoint was removed. See if face still exists.
                    label facei = reverseFaceMap[oldFacei];

                    if (facei >= 0)
                    {
                        if (isInternalFace(facei))
                        {
                            phi[facei] = phiU[facei];
                        }
                        else
                        {
                            label patchi = boundaryMesh().whichPatch(facei);
                            label i = facei - boundaryMesh()[patchi].start();

                            const fvsPatchScalarField& patchPhiU =
                                phiU.boundaryField()[patchi];
                            fvsPatchScalarField& patchPhi = phiBf[patchi];
                            patchPhi[i] = patchPhiU[i];
                        }
                    }
                }
            }
        }
    }


    // Update numbering of cells/vertices.
    meshCutter_.updateMesh(map);

    // Update numbering of protectedCells_
    if (protectedCells_.size())
    {
        PackedBoolList newProtectedCell(nCells());

        forAll(newProtectedCell, celli)
        {
            label oldCelli = map().cellMap()[celli];
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


const Foam::cellZone& Foam::dynamicRefineFvMesh::findCellZone
(
    const word& cellZoneName
) const
{
    const label cellZoneID = cellZones().findZoneID(cellZoneName);
    bool cellZoneFound = (cellZoneID != -1);
    reduce(cellZoneFound, orOp<bool>());
    if (!cellZoneFound)
    {
        FatalErrorInFunction
            << "cannot find cellZone " << cellZoneName
            << exit(FatalError);
    }

    return cellZones()[cellZoneID];
}


Foam::scalarField
Foam::dynamicRefineFvMesh::cellToPoint(const scalarField& vFld) const
{
    scalarField pFld(nPoints());

    forAll(pointCells(), pointi)
    {
        const labelList& pCells = pointCells()[pointi];

        scalar sum = 0.0;
        forAll(pCells, i)
        {
            sum += vFld[pCells[i]];
        }
        pFld[pointi] = sum/pCells.size();
    }
    return pFld;
}


Foam::scalarField Foam::dynamicRefineFvMesh::error
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


Foam::scalarField Foam::dynamicRefineFvMesh::error
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


void Foam::dynamicRefineFvMesh::selectRefineCandidates
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

    const labelList& cellLevel = meshCutter_.cellLevel();

    // Mark cells that are candidates for refinement.
    forAll(cellError, celli)
    {
        if
        (
            cellLevel[celli] < maxRefinement
         && cellError[celli] > 0
        )
        {
            candidateCells.set(celli, 1);
        }
    }
}


void Foam::dynamicRefineFvMesh::selectRefineCandidates
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

    const labelList& cellLevel = meshCutter_.cellLevel();

    // Mark cells that are candidates for refinement.
    forAll(cellError, celli)
    {
        if
        (
            cellLevel[celli] < maxRefinement
         && cellError[celli] > 0
        )
        {
            candidateCells.set(celli, 1);
        }
    }
}


Foam::scalar Foam::dynamicRefineFvMesh::selectRefineCandidates
(
    PackedBoolList& candidateCells,
    const dictionary& refineDict
) const
{
    const word fieldName(refineDict.lookup("field"));

    const volScalarField& vFld = lookupObject<volScalarField>(fieldName);

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


Foam::labelList Foam::dynamicRefineFvMesh::selectRefineCells
(
    const label maxCells,
    const label maxRefinement,
    const PackedBoolList& candidateCells
) const
{
    // Every refined cell causes 7 extra cells
    label nTotToRefine = (maxCells - globalData().nTotalCells()) / 7;

    const labelList& cellLevel = meshCutter_.cellLevel();

    // Mark cells that cannot be refined since they would trigger refinement
    // of protected cells (since 2:1 cascade)
    PackedBoolList unrefineableCells;
    calculateProtectedCells(unrefineableCells);

    // Count current selection
    label nLocalCandidates = count(candidateCells, 1);
    label nCandidates = returnReduce(nLocalCandidates, sumOp<label>());

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
        << " cells for refinement out of " << globalData().nTotalCells()
        << "." << endl;

    return consistentSet;
}


void Foam::dynamicRefineFvMesh::selectUnrefineCandidates
(
    boolList& unrefineCandidates,
    const volScalarField& vFld,
    const scalar unrefineLevel
) const
{
    forAll(pointCells(), pointi)
    {
        const labelList& pCells = pointCells()[pointi];

        scalar maxVal = -great;
        forAll(pCells, i)
        {
            maxVal = max(maxVal, vFld[pCells[i]]);
        }

        unrefineCandidates[pointi] =
            unrefineCandidates[pointi] && maxVal < unrefineLevel;
    }
}


void Foam::dynamicRefineFvMesh::selectUnrefineCandidates
(
    boolList& unrefineCandidates,
    const volScalarField& vFld,
    const cellZone& cZone,
    const scalar unrefineLevel
) const
{
    const Map<label>& zoneMap(cZone.lookupMap());

    forAll(pointCells(), pointi)
    {
        const labelList& pCells = pointCells()[pointi];

        scalar maxVal = -great;
        forAll(pCells, i)
        {
            if (zoneMap.found(pCells[i]))
            {
                maxVal = max(maxVal, vFld[pCells[i]]);
            }
        }

        unrefineCandidates[pointi] =
            unrefineCandidates[pointi] && maxVal < unrefineLevel;
    }
}


void Foam::dynamicRefineFvMesh::selectUnrefineCandidates
(
    boolList& unrefineCandidates,
    const dictionary& refineDict
) const
{
    if (refineDict.found("unrefineLevel"))
    {
        const word fieldName(refineDict.lookup("field"));
        const volScalarField& vFld
        (
            lookupObject<volScalarField>(fieldName)
        );

        const scalar unrefineLevel =
            refineDict.lookup<scalar>("unrefineLevel");

        if (refineDict.found("cellZone"))
        {
            selectUnrefineCandidates
            (
                unrefineCandidates,
                vFld,
                findCellZone(refineDict.lookup("cellZone")),
                unrefineLevel
            );
        }
        else
        {
            selectUnrefineCandidates
            (
                unrefineCandidates,
                vFld,
                unrefineLevel
            );
        }
    }
}


Foam::labelList Foam::dynamicRefineFvMesh::selectUnrefinePoints
(
    const PackedBoolList& markedCell,
    const boolList& unrefineCandidates
) const
{
    // All points that can be unrefined
    const labelList splitPoints(meshCutter_.getSplitPoints());

    DynamicList<label> newSplitPoints(splitPoints.size());

    forAll(splitPoints, i)
    {
        label pointi = splitPoints[i];

        if (unrefineCandidates[pointi])
        {
            // Check that all cells are not marked
            const labelList& pCells = pointCells()[pointi];

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


void Foam::dynamicRefineFvMesh::extendMarkedCells
(
    PackedBoolList& markedCell
) const
{
    // Mark faces using any marked cell
    boolList markedFace(nFaces(), false);

    forAll(markedCell, celli)
    {
        if (markedCell.get(celli))
        {
            const cell& cFaces = cells()[celli];

            forAll(cFaces, i)
            {
                markedFace[cFaces[i]] = true;
            }
        }
    }

    syncTools::syncFaceList(*this, markedFace, orEqOp<bool>());

    // Update cells using any markedFace
    for (label facei = 0; facei < nInternalFaces(); facei++)
    {
        if (markedFace[facei])
        {
            markedCell.set(faceOwner()[facei], 1);
            markedCell.set(faceNeighbour()[facei], 1);
        }
    }
    for (label facei = nInternalFaces(); facei < nFaces(); facei++)
    {
        if (markedFace[facei])
        {
            markedCell.set(faceOwner()[facei], 1);
        }
    }
}


void Foam::dynamicRefineFvMesh::checkEightAnchorPoints
(
    PackedBoolList& protectedCell,
    label& nProtected
) const
{
    const labelList& cellLevel = meshCutter_.cellLevel();
    const labelList& pointLevel = meshCutter_.pointLevel();

    labelList nAnchorPoints(nCells(), 0);

    forAll(pointLevel, pointi)
    {
        const labelList& pCells = pointCells(pointi);

        forAll(pCells, pCelli)
        {
            label celli = pCells[pCelli];

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

Foam::dynamicRefineFvMesh::dynamicRefineFvMesh(const IOobject& io)
:
    dynamicFvMesh(io),
    meshCutter_(*this),
    dumpLevel_(false),
    nRefinementIterations_(0),
    protectedCells_(nCells(), 0)
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

    labelList nAnchors(nCells(), 0);

    label nProtected = 0;

    forAll(pointCells(), pointi)
    {
        const labelList& pCells = pointCells()[pointi];

        forAll(pCells, i)
        {
            label celli = pCells[i];

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
        labelList neiLevel(nFaces());

        for (label facei = 0; facei < nInternalFaces(); facei++)
        {
            neiLevel[facei] = cellLevel[faceNeighbour()[facei]];
        }
        for (label facei = nInternalFaces(); facei < nFaces(); facei++)
        {
            neiLevel[facei] = cellLevel[faceOwner()[facei]];
        }
        syncTools::swapFaceList(*this, neiLevel);


        boolList protectedFace(nFaces(), false);

        forAll(faceOwner(), facei)
        {
            label faceLevel = max
            (
                cellLevel[faceOwner()[facei]],
                neiLevel[facei]
            );

            const face& f = faces()[facei];

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

        syncTools::syncFaceList(*this, protectedFace, orEqOp<bool>());

        for (label facei = 0; facei < nInternalFaces(); facei++)
        {
            if (protectedFace[facei])
            {
                protectedCells_.set(faceOwner()[facei], 1);
                nProtected++;
                protectedCells_.set(faceNeighbour()[facei], 1);
                nProtected++;
            }
        }
        for (label facei = nInternalFaces(); facei < nFaces(); facei++)
        {
            if (protectedFace[facei])
            {
                protectedCells_.set(faceOwner()[facei], 1);
                nProtected++;
            }
        }

        // Also protect any cells that are less than hex
        forAll(cells(), celli)
        {
            const cell& cFaces = cells()[celli];

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
                    if (faces()[cFaces[cFacei]].size() < 4)
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

        cellSet protectedCells(*this, "protectedCells", nProtected);
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

Foam::dynamicRefineFvMesh::~dynamicRefineFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dynamicRefineFvMesh::update()
{
    // Re-read dictionary. Chosen since usually -small so trivial amount
    // of time compared to actual refinement. Also very useful to be able
    // to modify on-the-fly.
    const dictionary refineDict
    (
        dynamicMeshDict().optionalSubDict(typeName + "Coeffs")
    );

    label refineInterval = refineDict.lookup<label>("refineInterval");

    bool hasChanged = false;

    if (refineInterval == 0)
    {
        topoChanging(hasChanged);

        return false;
    }
    else if (refineInterval < 0)
    {
        FatalErrorInFunction
            << "Illegal refineInterval " << refineInterval << nl
            << "The refineInterval setting in the dynamicMeshDict should"
            << " be >= 1." << nl
            << exit(FatalError);
    }

    // Note: cannot refine at time 0 since no V0 present since mesh not
    //       moved yet.

    if (time().timeIndex() > 0 && time().timeIndex() % refineInterval == 0)
    {
        label maxCells = refineDict.lookup<label>("maxCells");

        if (maxCells <= 0)
        {
            FatalErrorInFunction
                << "Illegal maximum number of cells " << maxCells << nl
                << "The maxCells setting in the dynamicMeshDict should"
                << " be > 0." << nl
                << exit(FatalError);
        }

        const label nBufferLayers =
            refineDict.lookup<label>("nBufferLayers");

        // Cells marked for refinement or otherwise protected from unrefinement.
        PackedBoolList refineCells(nCells());

        label maxRefinement = 0;

        if (refineDict.isDict("refinementRegions"))
        {
            const dictionary& refinementRegions
            (
                refineDict.subDict("refinementRegions")
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
            maxRefinement = selectRefineCandidates(refineCells, refineDict);
        }

        if (globalData().nTotalCells() < maxCells)
        {
            // Select subset of candidates. Take into account max allowable
            // cells, refinement level, protected cells.
            labelList cellsToRefine
            (
                selectRefineCells
                (
                    maxCells,
                    maxRefinement,
                    refineCells
                )
            );

            label nCellsToRefine = returnReduce
            (
                cellsToRefine.size(), sumOp<label>()
            );

            if (nCellsToRefine > 0)
            {
                // Refine/update mesh and map fields
                autoPtr<mapPolyMesh> map = refine(cellsToRefine);

                // Update refineCells. Note that some of the marked ones have
                // not been refined due to constraints.
                {
                    const labelList& cellMap = map().cellMap();
                    const labelList& reverseCellMap = map().reverseCellMap();

                    PackedBoolList newRefineCell(cellMap.size());

                    forAll(cellMap, celli)
                    {
                        label oldCelli = cellMap[celli];

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
                            newRefineCell.set(celli, refineCells.get(oldCelli));
                        }
                    }
                    refineCells.transfer(newRefineCell);
                }

                // Extend with a buffer layer to prevent neighbouring points
                // being unrefined.
                for (label i = 0; i < nBufferLayers; i++)
                {
                    extendMarkedCells(refineCells);
                }

                hasChanged = true;
            }
        }

        boolList unrefineCandidates(nPoints(), true);

        if (refineDict.isDict("refinementRegions"))
        {
            const dictionary& refinementRegions
            (
                refineDict.subDict("refinementRegions")
            );

            forAllConstIter(dictionary, refinementRegions, iter)
            {
                selectUnrefineCandidates
                (
                    unrefineCandidates,
                    refinementRegions.subDict(iter().keyword())
                );
            }
        }
        else
        {
            selectUnrefineCandidates
            (
                unrefineCandidates,
                refineDict
            );
        }

        {
            // Select unrefineable points that are not marked in refineCells
            labelList pointsToUnrefine
            (
                selectUnrefinePoints
                (
                    refineCells,
                    unrefineCandidates
                )
            );

            label nSplitPoints = returnReduce
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

    topoChanging(hasChanged);
    if (hasChanged)
    {
        // Reset moving flag (if any). If not using inflation we'll not move,
        // if are using inflation any follow on movePoints will set it.
        moving(false);
    }

    return hasChanged;
}


bool Foam::dynamicRefineFvMesh::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool write
) const
{
    // Force refinement data to go to the current time directory.
    const_cast<hexRef8&>(meshCutter_).setInstance(time().timeName());

    bool writeOk =
    (
        dynamicFvMesh::writeObject(fmt, ver, cmp, write)
     && meshCutter_.write(write)
    );

    if (dumpLevel_)
    {
        volScalarField scalarCellLevel
        (
            IOobject
            (
                "cellLevel",
                time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            *this,
            dimensionedScalar(dimless, 0)
        );

        const labelList& cellLevel = meshCutter_.cellLevel();

        forAll(cellLevel, celli)
        {
            scalarCellLevel[celli] = cellLevel[celli];
        }

        writeOk = writeOk && scalarCellLevel.write();
    }

    return writeOk;
}


// ************************************************************************* //
