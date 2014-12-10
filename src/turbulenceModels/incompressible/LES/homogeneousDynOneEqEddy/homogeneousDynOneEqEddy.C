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

#include "homogeneousDynOneEqEddy.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
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
    nuSgs_ = ck(D)*sqrt(k_)*delta();
    nuSgs_.correctBoundaryConditions();
}


dimensionedScalar homogeneousDynOneEqEddy::ck
(
    const volSymmTensorField& D
) const
{
    tmp<volScalarField> KK = 0.5*(filter_(magSqr(U())) - magSqr(filter_(U())));

    const volSymmTensorField MM
    (
        delta()*(filter_(sqrt(k_)*D) - 2*sqrt(KK + filter_(k_))*filter_(D))
    );

    dimensionedScalar MMMM = average(magSqr(MM));

    if (MMMM.value() > VSMALL)
    {
        tmp<volSymmTensorField> LL = dev(filter_(sqr(U())) - sqr(filter_(U())));

        return average(LL && MM)/MMMM;
    }
    else
    {
        return 0.0;
    }
}


dimensionedScalar homogeneousDynOneEqEddy::ce
(
    const volSymmTensorField& D
) const
{
    const volScalarField KK
    (
        0.5*(filter_(magSqr(U())) - magSqr(filter_(U())))
    );

    const volScalarField mm
    (
        pow(KK + filter_(k_), 1.5)/(2*delta()) - filter_(pow(k_, 1.5))/delta()
    );

    dimensionedScalar mmmm = average(magSqr(mm));

    if (mmmm.value() > VSMALL)
    {
        tmp<volScalarField> ee =
        (
            2*delta()*ck(D)
          * (
                filter_(sqrt(k_)*magSqr(D))
              - 2*sqrt(KK + filter_(k_))*magSqr(filter_(D))
            )
        );

        return average(ee*mm)/mmmm;
    }
    else
    {
        return 0.0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

homogeneousDynOneEqEddy::homogeneousDynOneEqEddy
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    LESModel(modelName, U, phi, transport, turbulenceModelName),
    GenEddyVisc(U, phi, transport),

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
    bound(k_, kMin_);

    updateSubGridScaleFields(symm(fvc::grad(U)));

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void homogeneousDynOneEqEddy::correct(const tmp<volTensorField>& tgradU)
{
    const volTensorField& gradU = tgradU();

    GenEddyVisc::correct(gradU);

    const volSymmTensorField D(symm(gradU));

    const volScalarField P(2.0*nuSgs_*magSqr(D));

    tmp<fvScalarMatrix> kEqn
    (
       fvm::ddt(k_)
     + fvm::div(phi(), k_)
     - fvm::laplacian(DkEff(), k_)
    ==
       P
     - fvm::Sp(ce(D)*sqrt(k_)/delta(), k_)
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
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
