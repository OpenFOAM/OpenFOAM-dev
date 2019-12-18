/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2019 OpenFOAM Foundation
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

#include "correctorConvergenceControl.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(correctorConvergenceControl, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::correctorConvergenceControl::getNSolves
(
    const fvMesh& mesh,
    const word& fieldName,
    label& n
)
{
    getNTypeSolves<scalar>(mesh, fieldName, n);
    getNTypeSolves<vector>(mesh, fieldName, n);
    getNTypeSolves<sphericalTensor>(mesh, fieldName, n);
    getNTypeSolves<symmTensor>(mesh, fieldName, n);
    getNTypeSolves<tensor>(mesh, fieldName, n);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::correctorConvergenceControl::correctorConvergenceControl
(
    const solutionControl& control,
    const word& loopName
)
:
    control_(control),
    loopName_(loopName)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::correctorConvergenceControl::~correctorConvergenceControl()
{}


// ************************************************************************* //
