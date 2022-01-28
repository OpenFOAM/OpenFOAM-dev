/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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

Class
    Foam::correctContactAngle

Description
    Correction for the boundary condition on the unit normal nHat on
    walls to produce the correct contact angle. The dynamic contact angle is
    calculated from the component of the velocity on the direction of the
    interface, parallel to the wall.

SourceFiles
    correctContactAngle.C

\*---------------------------------------------------------------------------*/

#ifndef correctContactAngle_H
#define correctContactAngle_H

#include "alphaContactAngleFvPatchScalarField.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Correct the contact angle for the two volume fraction fields
void correctContactAngle
(
    const volScalarField& alpha1,
    const volScalarField& alpha2,
    const volVectorField::Boundary& Ubf,
    const dimensionedScalar& deltaN,
    surfaceVectorField::Boundary& nHatbf
);

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
