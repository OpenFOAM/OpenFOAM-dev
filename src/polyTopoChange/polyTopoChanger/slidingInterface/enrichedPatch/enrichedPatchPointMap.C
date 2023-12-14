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

#include "enrichedPatch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::enrichedPatch::completePointMap() const
{
    if (pointMapComplete_)
    {
        FatalErrorInFunction
            << "Point map already completed"
            << abort(FatalError);
    }

    pointMapComplete_ = true;

    const Map<label>& pmm = pointMergeMap();

    // Get the mesh points for both patches.  If the point has not been
    // merged away, add it to the map

    // Do master patch
    const labelList& masterMeshPoints = masterPatch_.meshPoints();
    const pointField& masterLocalPoints = masterPatch_.localPoints();

    forAll(masterMeshPoints, pointi)
    {
        if (!pmm.found(masterMeshPoints[pointi]))
        {
            pointMap_.insert
            (
                masterMeshPoints[pointi],
                masterLocalPoints[pointi]
            );
        }
    }

    // Do slave patch
    const labelList& slaveMeshPoints = slavePatch_.meshPoints();
    const pointField& slaveLocalPoints = slavePatch_.localPoints();

    forAll(slaveMeshPoints, pointi)
    {
        if (!pmm.found(slaveMeshPoints[pointi]))
        {
            pointMap_.insert
            (
                slaveMeshPoints[pointi],
                slaveLocalPoints[pointi]
            );
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::Map<Foam::point>& Foam::enrichedPatch::pointMap()
{
    if (!pointMapComplete_)
    {
        completePointMap();
    }

    return pointMap_;
}


const Foam::Map<Foam::point>& Foam::enrichedPatch::pointMap() const
{
    if (!pointMapComplete_)
    {
        completePointMap();
    }

    return pointMap_;
}


// ************************************************************************* //
