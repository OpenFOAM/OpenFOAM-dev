/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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

#include "incompressibleTwoPhaseVoFMixture.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(incompressibleTwoPhaseVoFMixture, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::incompressibleTwoPhaseVoFMixture::incompressibleTwoPhaseVoFMixture
(
    const fvMesh& mesh
)
:
    twoPhaseVoFMixture(mesh),

    nuModel1_(viscosityModel::New(mesh, phase1Name())),
    nuModel2_(viscosityModel::New(mesh, phase2Name())),

    rho1_
    (
        IOobject
        (
            IOobject::groupName("rho", phase1Name()),
            mesh.time().constant(),
            mesh
        ),
        dimensionedScalar("rho", dimDensity, nuModel1_())
    ),
    rho2_
    (
        IOobject
        (
            IOobject::groupName("rho", phase2Name()),
            mesh.time().constant(),
            mesh
        ),
        dimensionedScalar("rho", dimDensity, nuModel2_())
    ),

    rho_
    (
        IOobject
        (
            "rho",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("rho", dimDensity, 0)
    ),

    nu_
    (
        IOobject
        (
            "nu",
            mesh.time().name(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimKinematicViscosity, 0),
        calculatedFvPatchScalarField::typeName
    )
{
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::incompressibleTwoPhaseVoFMixture::correct()
{
    rho_ = alpha1()*rho1_ + alpha2()*rho2_;

    nuModel1_->correct();
    nuModel2_->correct();

    const volScalarField limitedAlpha1
    (
        "limitedAlpha1",
        min(max(alpha1(), scalar(0)), scalar(1))
    );

    // Average kinematic viscosity calculated from dynamic viscosity
    nu_ =
    (
        limitedAlpha1*rho1_*nuModel1_->nu()
      + (scalar(1) - limitedAlpha1)*rho2_*nuModel2_->nu()
    )/(limitedAlpha1*rho1_ + (scalar(1) - limitedAlpha1)*rho2_);
}


bool Foam::incompressibleTwoPhaseVoFMixture::read()
{
    if (twoPhaseVoFMixture::read())
    {
        if (nuModel1_->read() && nuModel2_->read())
        {
            nuModel1_->lookup("rho") >> rho1_;
            nuModel2_->lookup("rho") >> rho2_;

            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
