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

\*---------------------------------------------------------------------------*/

#include "cavitatingTwoPhaseMixture.H"
#include "barotropicCompressibilityModel.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cavitatingTwoPhaseMixture, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cavitatingTwoPhaseMixture::cavitatingTwoPhaseMixture
(
    const fvMesh& mesh
)
:
    twoPhaseMixture(mesh),

    alphav_(alpha1()),
    alphal_(alpha2()),

    nuModelv_(viscosityModel::New(mesh, phase1Name())),
    nuModell_(viscosityModel::New(mesh, phase2Name())),

    rhov_("rho", dimDensity, nuModelv_()),
    rhol_("rho", dimDensity, nuModell_()),

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
    ),

    thermodynamicProperties_
    (
        IOobject
        (
            "thermodynamicProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    psil_
    (
        "psil",
        dimCompressibility,
        thermodynamicProperties_
    ),

    rholSat_
    (
        "rholSat",
        dimDensity,
        thermodynamicProperties_
    ),

    psiv_
    (
        "psiv",
        dimCompressibility,
        thermodynamicProperties_
    ),

    pSat_
    (
        "pSat",
        dimPressure,
        thermodynamicProperties_
    ),

    rhovSat_("rhovSat", psiv_*pSat_),

    rhol0_("rhol0", rholSat_ - pSat_*psil_),

    rhoMin_
    (
        "rhoMin",
        dimDensity,
        thermodynamicProperties_
    ),

    p_
    (
        IOobject
        (
            "p",
            mesh.time().name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    psiModel_
    (
        barotropicCompressibilityModel::New
        (
            thermodynamicProperties_,
            alphav_
        )
    ),

    psi_(psiModel_->psi()),

    rho_
    (
        IOobject
        (
            "rho",
            mesh.time().name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    )
{
    alphav_.oldTime();

    rho_ = max
    (
        psi_*p_
      + alphal_*rhol0_
      + ((alphav_*psiv_ + alphal_*psil_) - psi_)*pSat_,
        rhoMin_
    );

    correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cavitatingTwoPhaseMixture::~cavitatingTwoPhaseMixture()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::cavitatingTwoPhaseMixture::read()
{
    if (twoPhaseMixture::read())
    {
        nuModelv_->lookup("rho") >> rhov_;
        nuModell_->lookup("rho") >> rhol_;

        return true;
    }
    else
    {
        return false;
    }
}


void Foam::cavitatingTwoPhaseMixture::correctPressure()
{
    p_ =
    (
        rho_
      - alphal_*rhol0_
      - ((alphav_*psiv_ + alphal_*psil_) - psi_)*pSat_
    )/psi_;

    p_.correctBoundaryConditions();
}


void Foam::cavitatingTwoPhaseMixture::correct()
{
    rho_ == max(rho_, rhoMin_);

    alphav_ =
        max
        (
            min
            (
                (rho_ - rholSat_)/(rhovSat_ - rholSat_),
                scalar(1)
            ),
            scalar(0)
        );
    alphal_ = 1.0 - alphav_;

    Info<< "max-min alphav: " << max(alphav_).value()
        << " " << min(alphav_).value() << endl;

    psiModel_->correct();

    nuModelv_->correct();
    nuModell_->correct();

    const volScalarField limitedAlphav
    (
        "limitedAlphav",
        min(max(alphav_, scalar(0)), scalar(1))
    );

    // Mixture kinematic viscosity calculated from mixture dynamic viscosity
    nu_ =
    (
        limitedAlphav*rhov_*nuModelv_->nu()
      + (scalar(1) - limitedAlphav)*rhol_*nuModell_->nu()
    )/(limitedAlphav*rhov_ + (scalar(1) - limitedAlphav)*rhol_);
}


// ************************************************************************* //
