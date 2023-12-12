/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "attachDetach.H"
#include "polyMesh.H"
#include "primitiveMesh.H"
#include "primitiveFacePatch.H"
#include "polyTopoChanger.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::Map<Foam::label>&
Foam::attachDetach::pointMatchMap() const
{
    if (!pointMatchMapPtr_)
    {
        calcPointMatchMap();
    }

    return *pointMatchMapPtr_;
}


void Foam::attachDetach::calcPointMatchMap() const
{
    if (debug)
    {
        Pout<< "void attachDetach::calcPointMatchMap() const "
            << " for object " << name() << " : "
            << "Calculating point matching" << endl;
    }

    if (pointMatchMapPtr_)
    {
        FatalErrorInFunction
            << "Point match map already calculated for object " << name()
            << abort(FatalError);
    }

    const polyMesh& mesh = topoChanger().mesh();
    const faceList& faces = mesh.faces();

    const polyPatch& masterPatch = mesh.boundaryMesh()[masterPatchID_.index()];
    const polyPatch& slavePatch = mesh.boundaryMesh()[slavePatchID_.index()];

    // Create the reverse patch out of the slave patch
    primitiveFacePatch reverseSlavePatch
    (
        faceList(slavePatch.size()),
        mesh.points()
    );

    const label slavePatchStart = slavePatch.start();

    forAll(reverseSlavePatch, facei)
    {
        reverseSlavePatch[facei] =
            faces[slavePatchStart + facei].reverseFace();
    }

    // Create point merge list and remove merged points
    const labelList& masterMeshPoints = masterPatch.meshPoints();
    const labelList& slaveMeshPoints = reverseSlavePatch.meshPoints();

    const faceList& masterLocalFaces = masterPatch.localFaces();
    const faceList& slaveLocalFaces = reverseSlavePatch.localFaces();

    pointMatchMapPtr_ = new Map<label>(2*slaveMeshPoints.size());
    Map<label>& removedPointMap = *pointMatchMapPtr_;

    forAll(masterLocalFaces, facei)
    {
        const face& curMasterPoints = masterLocalFaces[facei];
        const face& curSlavePoints = slaveLocalFaces[facei];

        forAll(curMasterPoints, pointi)
        {
            // If the master and slave point labels are the same, the
            // point remains.  Otherwise, the slave point is removed and
            // replaced by the master
            if
            (
                masterMeshPoints[curMasterPoints[pointi]]
             != slaveMeshPoints[curSlavePoints[pointi]]
            )
            {
                // Pout<< "Matching slave point "
                //     << slaveMeshPoints[curSlavePoints[pointi]]
                //     << " with "
                //     << masterMeshPoints[curMasterPoints[pointi]]
                //     << endl;

                // Grab the addressing
                removedPointMap.insert
                (
                    slaveMeshPoints[curSlavePoints[pointi]],
                    masterMeshPoints[curMasterPoints[pointi]]
                );
            }
        }
    }

    if (debug)
    {
        Pout<< "void attachDetach::calcPointMatchMap() const "
            << " for object " << name() << " : "
            << "Finished calculating point matching" << endl;
    }
}


// ************************************************************************* //
