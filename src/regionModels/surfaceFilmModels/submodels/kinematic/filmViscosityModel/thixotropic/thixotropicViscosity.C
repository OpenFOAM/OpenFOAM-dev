/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2022 OpenFOAM Foundation
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

#include "thixotropicViscosity.H"
#include "kinematicSingleLayer.H"
#include "addToRunTimeSelectionTable.H"

#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvcDiv.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(thixotropicViscosity, 0);

addToRunTimeSelectionTable
(
    viscosityModel,
    thixotropicViscosity,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

thixotropicViscosity::thixotropicViscosity
(
    surfaceFilmRegionModel& film,
    const dictionary& dict,
    volScalarField& mu
)
:
    viscosityModel(typeName, film, dict, mu),
    a_("a", dimless/dimTime, coeffDict_),
    b_("b", dimless, coeffDict_),
    d_("d", dimless, coeffDict_),
    c_("c", pow(dimTime, d_.value() - scalar(1)), coeffDict_),
    mu0_("mu0", dimPressure*dimTime, coeffDict_),
    muInf_("muInf", mu0_.dimensions(), coeffDict_),
    BinghamPlastic_(coeffDict_.found("tauy")),
    tauy_
    (
        BinghamPlastic_
      ? dimensionedScalar("tauy", dimPressure, coeffDict_)
      : dimensionedScalar("tauy", dimPressure, 0)
    ),
    K_(1 - sqrt(muInf_/mu0_)),
    lambda_
    (
        IOobject
        (
            typedName("lambda"),
            film.regionMesh().time().timeName(),
            film.regionMesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        film.regionMesh()
    )
{
    lambda_.min(1);
    lambda_.max(0);

    // Initialise viscosity to inf value because it cannot be evaluated yet
    mu_ = muInf_;
    mu_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

thixotropicViscosity::~thixotropicViscosity()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void thixotropicViscosity::correct
(
    const volScalarField& p,
    const volScalarField& T
)
{
    const kinematicSingleLayer& film = filmType<kinematicSingleLayer>();

    const volScalarField::Internal& delta = film.delta();
    const volScalarField::Internal alphaRho(film.alpha()()*film.rho()());
    const volScalarField::Internal& coverage = film.coverage();
    const surfaceScalarField& phiU = film.phiU();

    // Shear rate
    const volScalarField::Internal gDot
    (
        "gDot",
        coverage*mag(film.Us() - film.Uw())/(delta + film.deltaSmall())
    );

    const dimensionedScalar alphaRho0
    (
        "alphaRho0",
        alphaRho.dimensions(),
        rootVSmall
    );

    fvScalarMatrix lambdaEqn
    (
        fvm::ddt(lambda_)
      + fvm::div(phiU, lambda_)
      - fvm::Sp(fvc::div(phiU), lambda_)
      ==
        a_*pow((1 - lambda_), b_)
      - fvm::Sp(c_*pow(gDot, d_), lambda_)

        // Include the effect of the impinging droplets added with lambda = 0
      - fvm::Sp
        (
            max
            (
               -film.rhoSp(),
                dimensionedScalar(film.rhoSp().dimensions(), 0)
            )/(alphaRho + alphaRho0),
            lambda_
        )
    );

    lambdaEqn.relax();
    lambdaEqn.solve();

    lambda_.min(1);
    lambda_.max(0);

    mu_ = muInf_/(sqr(1 - K_*lambda_) + rootVSmall);

    // Add optional yield stress contribution to the viscosity
    if (BinghamPlastic_)
    {
        dimensionedScalar tauySmall("tauySmall", tauy_.dimensions(), small);
        dimensionedScalar muMax_("muMax", 100*mu0_);

        mu_.ref() = min
        (
            tauy_/(gDot + 1.0e-4*(tauy_ + tauySmall)/mu0_) + mu_(),
            muMax_
        );
    }

    mu_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
