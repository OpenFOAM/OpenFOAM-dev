/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2024 OpenFOAM Foundation
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

#include "preserveFaceZonesConstraint.H"
#include "addToRunTimeSelectionTable.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace decompositionConstraints
{
    defineTypeName(preserveFaceZonesConstraint);

    addToRunTimeSelectionTable
    (
        decompositionConstraint,
        preserveFaceZonesConstraint,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::decompositionConstraints::preserveFaceZonesConstraint::
preserveFaceZonesConstraint
(
    const dictionary& constraintsDict,
    const word& modelType
)
:
    decompositionConstraint(constraintsDict, typeName),
    zones_(constraintsDict.lookup("zones"))
{
    if (decompositionConstraint::debug)
    {
        Info<< type() << " : adding constraints to keep owner and neighbour"
            << " of faces in zones " << zones_
            << " on same processor" << endl;
    }
}


Foam::decompositionConstraints::preserveFaceZonesConstraint::
preserveFaceZonesConstraint
(
    const wordReList& zones
)
:
    decompositionConstraint(dictionary(), typeName),
    zones_(zones)
{
    if (decompositionConstraint::debug)
    {
        Info<< type() << " : adding constraints to keep owner and neighbour"
            << " of faces in zones " << zones_
            << " on same processor" << endl;
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::decompositionConstraints::preserveFaceZonesConstraint::add
(
    const polyMesh& mesh,
    boolList& blockedFace,
    PtrList<labelList>& specifiedProcessorFaces,
    labelList& specifiedProcessor,
    List<labelPair>& explicitConnections
) const
{
    blockedFace.setSize(mesh.nFaces(), true);

    const faceZoneList& fZones = mesh.faceZones();

    const labelList zoneIDs = findStrings(zones_, fZones.toc());

    label nUnblocked = 0;

    forAll(zoneIDs, i)
    {
        const faceZone& fz = fZones[zoneIDs[i]];

        forAll(fz, i)
        {
            if (blockedFace[fz[i]])
            {
                blockedFace[fz[i]] = false;
                nUnblocked++;
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


void Foam::decompositionConstraints::preserveFaceZonesConstraint::apply
(
    const polyMesh& mesh,
    const boolList& blockedFace,
    const PtrList<labelList>& specifiedProcessorFaces,
    const labelList& specifiedProcessor,
    const List<labelPair>& explicitConnections,
    labelList& decomposition
) const
{
    // If the decomposition has not enforced the constraint do it over
    // here.


    // Synchronise decomposition on boundary
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    labelList destProc(mesh.nFaces()-mesh.nInternalFaces(), labelMax);

    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];

        const labelUList& faceCells = pp.faceCells();

        forAll(faceCells, i)
        {
            label bFaceI = pp.start()+i-mesh.nInternalFaces();
            destProc[bFaceI] = decomposition[faceCells[i]];
        }
    }

    syncTools::syncBoundaryFaceList(mesh, destProc, minEqOp<label>());


    // Override if differing
    // ~~~~~~~~~~~~~~~~~~~~~

    const faceZoneList& fZones = mesh.faceZones();

    const labelList zoneIDs = findStrings(zones_, fZones.toc());

    label nChanged = 0;

    forAll(zoneIDs, i)
    {
        const faceZone& fz = fZones[zoneIDs[i]];

        forAll(fz, i)
        {
            label faceI = fz[i];

            label own = mesh.faceOwner()[faceI];

            if (mesh.isInternalFace(faceI))
            {
                label nei = mesh.faceNeighbour()[faceI];
                if (decomposition[own] != decomposition[nei])
                {
                    decomposition[nei] = decomposition[own];
                    nChanged++;
                }
            }
            else
            {
                label bFaceI = faceI-mesh.nInternalFaces();
                if (decomposition[own] != destProc[bFaceI])
                {
                    decomposition[own] = destProc[bFaceI];
                    nChanged++;
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
