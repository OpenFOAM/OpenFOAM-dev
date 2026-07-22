/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2026 OpenFOAM Foundation
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

#include "singleProcessorFaceZonesConstraint.H"
#include "addToRunTimeSelectionTable.H"
#include "syncTools.H"
#include "faceSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace decompositionConstraints
{
    defineTypeName(singleProcessorFaceZonesConstraint);

    addToRunTimeSelectionTable
    (
        decompositionConstraint,
        singleProcessorFaceZonesConstraint,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decompositionConstraints::singleProcessorFaceZonesConstraint::
singleProcessorFaceZonesConstraint
(
    const dictionary& constraintsDict,
    const word& modelType
)
:
    decompositionConstraint(constraintsDict, typeName),
    zoneNameAndProcs_(constraintsDict.lookup("singleProcessorFaceZones"))
{
    if (decompositionConstraint::debug)
    {
        Info<< type()
            << " : adding constraints to keep" << endl;

        forAll(zoneNameAndProcs_, zonei)
        {
            Info<< "    all cells connected to faceZone "
                << zoneNameAndProcs_[zonei].first() << " on processor "
                << zoneNameAndProcs_[zonei].second() << endl;
        }
    }
}


Foam::decompositionConstraints::singleProcessorFaceZonesConstraint::
singleProcessorFaceZonesConstraint
(
    const List<Tuple2<word, label>>& zoneNameAndProcs
)
:
    decompositionConstraint(dictionary(), typeName),
    zoneNameAndProcs_(zoneNameAndProcs)
{
    if (decompositionConstraint::debug)
    {
        Info<< type()
            << " : adding constraints to keep" << endl;

        forAll(zoneNameAndProcs_, zonei)
        {
            Info<< "    all cells connected to faceZone "
                << zoneNameAndProcs_[zonei].first() << " on processor "
                << zoneNameAndProcs_[zonei].second() << endl;
        }
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::decompositionConstraints::singleProcessorFaceZonesConstraint::add
(
    const polyMesh& mesh,
    boolList& blockedFace,
    PtrList<labelList>& specifiedProcessorFaces,
    labelList& specifiedProcessor,
    List<labelPair>& explicitConnections
) const
{
    blockedFace.setSize(mesh.nFaces(), true);

    // Mark faces already in set
    labelList faceToSet(mesh.nFaces(), -1);
    forAll(specifiedProcessorFaces, zonei)
    {
        const labelList& faceLabels = specifiedProcessorFaces[zonei];
        forAll(faceLabels, i)
        {
            faceToSet[faceLabels[i]] = zonei;
        }
    }

    forAll(zoneNameAndProcs_, zonei)
    {
        const faceZone& fz =
            mesh.faceZones()[zoneNameAndProcs_[zonei].first()];

        // Check that it does not overlap with existing specifiedProcessorFaces
        labelList nMatch(specifiedProcessorFaces.size(), 0);
        forAll(fz, fzFacei)
        {
            const label zonei = faceToSet[fz[fzFacei]];
            if (zonei != -1)
            {
                nMatch[zonei]++;
            }
        }

        // Only store if all faces are not yet in specifiedProcessorFaces
        // (on all processors)
        bool store = true;
        forAll(nMatch, zonei)
        {
            if (nMatch[zonei] == fz.size())
            {
                // full match
                store = false;
                break;
            }
            else if (nMatch[zonei] > 0)
            {
                // partial match
                store = false;
                break;
            }
        }
        reduce(store, andOp<bool>());

        if (store)
        {
            specifiedProcessorFaces.append(fz.labelList::clone());
            specifiedProcessor.append(zoneNameAndProcs_[zonei].second());
        }
    }

    // Unblock all point connected faces
    // 1. Mark all points on specifiedProcessorFaces
    boolList procFacePoint(mesh.nPoints(), false);
    forAll(specifiedProcessorFaces, zonei)
    {
        const labelList& set = specifiedProcessorFaces[zonei];
        forAll(set, fI)
        {
            const face& f = mesh.faces()[set[fI]];
            forAll(f, fp)
            {
                procFacePoint[f[fp]] = true;
            }
        }
    }
    syncTools::syncPointList(mesh, procFacePoint, orEqOp<bool>(), false);

    // 2. Unblock all faces on procFacePoint
    label nUnblocked = 0;
    forAll(procFacePoint, pointi)
    {
        if (procFacePoint[pointi])
        {
            const labelList& pFaces = mesh.pointFaces()[pointi];
            forAll(pFaces, i)
            {
                if (blockedFace[pFaces[i]])
                {
                    blockedFace[pFaces[i]] = false;
                    nUnblocked++;
                }
            }
        }
    }

    if (decompositionConstraint::debug & 2)
    {
        reduce(nUnblocked, sumOp<label>());
        Info<< type() << " : unblocked " << nUnblocked << " faces" << endl;
    }

    syncTools::syncFaceList(mesh, blockedFace, andEqOp<bool>());
}


void Foam::decompositionConstraints::singleProcessorFaceZonesConstraint::apply
(
    const polyMesh& mesh,
    const boolList& blockedFace,
    const PtrList<labelList>& specifiedProcessorFaces,
    const labelList& specifiedProcessor,
    const List<labelPair>& explicitConnections,
    labelList& decomposition
) const
{
    // For specifiedProcessorFaces rework the cellToProc to enforce
    // all on one processor since we can't guarantee that the input
    // to regionSplit was a single region.
    // E.g. faceSet 'a' with the cells split into two regions
    // by a notch formed by two walls
    //
    //          \   /
    //           \ /
    //    ---a----+-----a-----
    //
    //
    // Note that reworking the cellToProc might make the decomposition
    // unbalanced.
    label nChanged = 0;

    forAll(specifiedProcessorFaces, zonei)
    {
        const labelList& set = specifiedProcessorFaces[zonei];

        // Get the processor to use for the set
        label procI = specifiedProcessor[zonei];
        if (procI == -1)
        {
            // If no processor specified use the one from the
            // 0th element
            if (set.size())
            {
                procI = decomposition[mesh.faceOwner()[set[0]]];
            }
            reduce(procI, maxOp<label>());
        }

        // Get all points on the sets
        boolList procFacePoint(mesh.nPoints(), false);
        forAll(set, fI)
        {
            const face& f = mesh.faces()[set[fI]];
            forAll(f, fp)
            {
                procFacePoint[f[fp]] = true;
            }
        }
        syncTools::syncPointList(mesh, procFacePoint, orEqOp<bool>(), false);

        // 2. Unblock all faces on procFacePoint
        forAll(procFacePoint, pointi)
        {
            if (procFacePoint[pointi])
            {
                const labelList& pFaces = mesh.pointFaces()[pointi];
                forAll(pFaces, i)
                {
                    label faceI = pFaces[i];

                    label own = mesh.faceOwner()[faceI];
                    if (decomposition[own] != procI)
                    {
                        decomposition[own] = procI;
                        nChanged++;
                    }
                    if (mesh.isInternalFace(faceI))
                    {
                        label nei = mesh.faceNeighbour()[faceI];
                        if (decomposition[nei] != procI)
                        {
                            decomposition[nei] = procI;
                            nChanged++;
                        }
                    }
                }
            }
        }
    }

    if (decompositionConstraint::debug & 2)
    {
        reduce(nChanged, sumOp<label>());
        Info<< type() << " : changed decomposition on " << nChanged
            << " cells" << endl;
    }
}


// ************************************************************************* //
