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

#include "incompressibleTwoPhaseVoFMixture.H"
#include "surfaceInterpolate.H"
#include "addToRunTimeSelectionTable.H"

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

    rho1_("rho", dimDensity, nuModel1_()),
    rho2_("rho", dimDensity, nuModel2_()),

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
        dimensionedScalar(dimViscosity, 0),
        calculatedFvPatchScalarField::typeName
    )
{
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::incompressibleTwoPhaseVoFMixture::mu() const
{
    const volScalarField limitedAlpha1
    (
        min(max(alpha1(), scalar(0)), scalar(1))
    );

    return volScalarField::New
    (
        "mu",
        limitedAlpha1*rho1_*nuModel1_->nu()
      + (scalar(1) - limitedAlpha1)*rho2_*nuModel2_->nu()
    );
}


Foam::tmp<Foam::surfaceScalarField>
Foam::incompressibleTwoPhaseVoFMixture::muf() const
{
    const surfaceScalarField alpha1f
    (
        min(max(fvc::interpolate(alpha1()), scalar(0)), scalar(1))
    );

    return surfaceScalarField::New
    (
        "muf",
        alpha1f*rho1_*fvc::interpolate(nuModel1_->nu())
      + (scalar(1) - alpha1f)*rho2_*fvc::interpolate(nuModel2_->nu())
    );
}


Foam::tmp<Foam::surfaceScalarField>
Foam::incompressibleTwoPhaseVoFMixture::nuf() const
{
    const surfaceScalarField alpha1f
    (
        min(max(fvc::interpolate(alpha1()), scalar(0)), scalar(1))
    );

    return surfaceScalarField::New
    (
        "nuf",
        (
            alpha1f*rho1_*fvc::interpolate(nuModel1_->nu())
          + (scalar(1) - alpha1f)*rho2_*fvc::interpolate(nuModel2_->nu())
        )/(alpha1f*rho1_ + (scalar(1) - alpha1f)*rho2_)
    );
}


bool Foam::incompressibleTwoPhaseVoFMixture::read()
{
    if (twoPhaseVoFMixture::read())
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
    nu_ = mu()/(limitedAlpha1*rho1_ + (scalar(1) - limitedAlpha1)*rho2_);
}


// ************************************************************************* //
