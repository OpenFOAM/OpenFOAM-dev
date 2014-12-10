/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "homogeneousDynOneEqEddy.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(homogeneousDynOneEqEddy, 0);
addToRunTimeSelectionTable(LESModel, homogeneousDynOneEqEddy, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void homogeneousDynOneEqEddy::updateSubGridScaleFields
(
    const volSymmTensorField& D
)
{
    muSgs_ = ck_(D)*rho()*sqrt(k_)*delta();
    muSgs_.correctBoundaryConditions();

    alphaSgs_ = muSgs_/Prt_;
    alphaSgs_.correctBoundaryConditions();
}


dimensionedScalar homogeneousDynOneEqEddy::ck_
(
    const volSymmTensorField& D
) const
{
    volScalarField KK(0.5*(filter_(magSqr(U())) - magSqr(filter_(U()))));

    volSymmTensorField LL(dev(filter_(sqr(U())) - (sqr(filter_(U())))));

    volSymmTensorField MM
    (
        delta()*(filter_(sqrt(k_)*D) - 2*sqrt(KK + filter_(k_))*filter_(D))
    );

    return average(LL && MM)/average(magSqr(MM));
}


dimensionedScalar homogeneousDynOneEqEddy::ce_
(
    const volSymmTensorField& D
) const
{
    volScalarField KK(0.5*(filter_(magSqr(U())) - magSqr(filter_(U()))));

    volScalarField mm
    (
        pow(KK + filter_(k_), 1.5)/(2*delta()) - filter_(pow(k_, 1.5))/delta()
    );

    volScalarField ee
    (
        2*delta()*ck_(D)
       *(
            filter_(sqrt(k_)*magSqr(D))
          - 2*sqrt(KK + filter_(k_))*magSqr(filter_(D))
        )
    );

    return average(ee*mm)/average(mm*mm);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

homogeneousDynOneEqEddy::homogeneousDynOneEqEddy
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const fluidThermo& thermoPhysicalModel,
    const word& turbulenceModelName,
    const word& modelName
)
:
    LESModel(modelName, rho, U, phi, thermoPhysicalModel, turbulenceModelName),
    GenEddyVisc(rho, U, phi, thermoPhysicalModel),

    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    filterPtr_(LESfilter::New(U.mesh(), coeffDict())),
    filter_(filterPtr_())
{
    updateSubGridScaleFields(dev(symm(fvc::grad(U))));

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void homogeneousDynOneEqEddy::correct(const tmp<volTensorField>& tgradU)
{
    const volTensorField& gradU = tgradU();

    GenEddyVisc::correct(gradU);

    volSymmTensorField D(dev(symm(gradU)));
    volScalarField divU(fvc::div(phi()/fvc::interpolate(rho())));
    volScalarField G(GName(), 2*muSgs_*(gradU && D));

    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(rho(), k_)
      + fvm::div(phi(), k_)
      - fvm::laplacian(DkEff(), k_)
     ==
        G
      - fvm::SuSp(2.0/3.0*rho()*divU, k_)
      - fvm::Sp(ce_(D)*rho()*sqrt(k_)/delta(), k_)
    );

    kEqn().relax();
    kEqn().solve();

    bound(k_, kMin_);

    updateSubGridScaleFields(D);
}


bool homogeneousDynOneEqEddy::read()
{
    if (GenEddyVisc::read())
    {
        filter_.read(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
