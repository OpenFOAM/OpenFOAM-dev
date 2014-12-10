/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "slidingInterface.H"
#include "polyTopoChange.H"
#include "polyMesh.H"
#include "polyTopoChanger.H"
#include "polyRemovePoint.H"
#include "polyRemoveFace.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::slidingInterface::clearCouple
(
    polyTopoChange& ref
) const
{
    if (debug)
    {
        Pout<< "void slidingInterface::clearCouple("
            << "polyTopoChange& ref) const for object " << name() << " : "
            << "Clearing old couple points and faces." << endl;
    }

    // Remove all points from the point zone

    const polyMesh& mesh = topoChanger().mesh();

    const labelList& cutPointZoneLabels =
        mesh.pointZones()[cutPointZoneID_.index()];

    forAll(cutPointZoneLabels, pointI)
    {
        ref.setAction(polyRemovePoint(cutPointZoneLabels[pointI]));
    }

    // Remove all faces from the face zone
    const labelList& cutFaceZoneLabels =
        mesh.faceZones()[cutFaceZoneID_.index()];

    forAll(cutFaceZoneLabels, faceI)
    {
        ref.setAction(polyRemoveFace(cutFaceZoneLabels[faceI]));
    }

    if (debug)
    {
        Pout<< "void slidingInterface::clearCouple("
            << "polyTopoChange& ref) const for object " << name() << " : "
            << "Finished clearing old couple points and faces." << endl;
    }
}


// ************************************************************************* //
