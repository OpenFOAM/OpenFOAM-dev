/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "motionDirectionalDiffusivity.H"
#include "surfaceInterpolate.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(motionDirectionalDiffusivity, 0);

    addToRunTimeSelectionTable
    (
        motionDiffusivity,
        motionDirectionalDiffusivity,
        Istream
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::motionDirectionalDiffusivity::motionDirectionalDiffusivity
(
    const fvMesh& mesh,
    Istream& mdData
)
:
    motionDiffusivity(mesh),
    diffusivityVector_(mdData)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::motionDirectionalDiffusivity::~motionDirectionalDiffusivity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField>
Foam::motionDirectionalDiffusivity::operator()() const
{
    if (mesh().foundObject<volVectorField>("cellMotionU"))
    {
        const volVectorField& cellMotionU =
            mesh().lookupObject<volVectorField>("cellMotionU");

        volVectorField D
        (
            IOobject
            (
                "D",
                mesh().time().name(),
                mesh()
            ),
            diffusivityVector_.y()*vector::one
          + (diffusivityVector_.x() - diffusivityVector_.y())*cellMotionU
           /(mag(cellMotionU) + dimensionedScalar(dimVelocity, small)),
            zeroGradientFvPatchVectorField::typeName
        );
        D.correctBoundaryConditions();

        const surfaceVectorField n(mesh().Sf()/mesh().magSf());

        return surfaceScalarField::New
        (
            "faceDiffusivity",
            (n & cmptMultiply(fvc::interpolate(D), n))
        );
    }
    else
    {
        return tmp<surfaceScalarField>
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "faceDiffusivity",
                    mesh().time().name(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh(),
                dimensionedScalar(dimless, 1.0)
            )
        );
    }
}


// ************************************************************************* //
