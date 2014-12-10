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

#include "cyclicSlipPointPatch.H"
#include "pointConstraint.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(cyclicSlipPointPatch, 0);

// Add the patch constructor functions to the hash tables
addToRunTimeSelectionTable
(
    facePointPatch,
    cyclicSlipPointPatch,
    polyPatch
);


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const vectorField& cyclicSlipPointPatch::pointNormals() const
{
    // Use underlying patch normals
    return refCast<const facePointPatch>
    (
        *this
    ).facePointPatch::pointNormals();
}


void cyclicSlipPointPatch::applyConstraint
(
    const label pointi,
    pointConstraint& pc
) const
{
    pc.applyConstraint(pointNormals()[pointi]);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
