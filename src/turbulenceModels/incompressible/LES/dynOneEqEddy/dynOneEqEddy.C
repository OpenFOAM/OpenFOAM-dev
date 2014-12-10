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

#include "dynOneEqEddy.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(dynOneEqEddy, 0);
addToRunTimeSelectionTable(LESModel, dynOneEqEddy, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void dynOneEqEddy::updateSubGridScaleFields
(
    const volSymmTensorField& D,
    const volScalarField& KK
)
{
    nuSgs_ = ck(D, KK)*sqrt(k_)*delta();
    nuSgs_.correctBoundaryConditions();
}


volScalarField dynOneEqEddy::ck
(
    const volSymmTensorField& D,
    const volScalarField& KK
) const
{
    const volSymmTensorField LL
    (
        simpleFilter_(dev(filter_(sqr(U())) - (sqr(filter_(U())))))
    );

    const volSymmTensorField MM
    (
        simpleFilter_(-2.0*delta()*sqrt(KK)*filter_(D))
    );

    const volScalarField ck
    (
        simpleFilter_(0.5*(LL && MM))
       /(
            simpleFilter_(magSqr(MM))
          + dimensionedScalar("small", sqr(MM.dimensions()), VSMALL)
        )
    );

    tmp<volScalarField> tfld = 0.5*(mag(ck) + ck);
    return tfld();
}


volScalarField dynOneEqEddy::ce
(
    const volSymmTensorField& D,
    const volScalarField& KK
) const
{
    const volScalarField ce
    (
        simpleFilter_(nuEff()*(filter_(magSqr(D)) - magSqr(filter_(D))))
       /simpleFilter_(pow(KK, 1.5)/(2.0*delta()))
    );

    tmp<volScalarField> tfld = 0.5*(mag(ce) + ce);
    return tfld();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

dynOneEqEddy::dynOneEqEddy
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

    simpleFilter_(U.mesh()),
    filterPtr_(LESfilter::New(U.mesh(), coeffDict())),
    filter_(filterPtr_())
{
    bound(k_, kMin_);

    const volScalarField KK(0.5*(filter_(magSqr(U)) - magSqr(filter_(U))));
    updateSubGridScaleFields(symm(fvc::grad(U)), KK);

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void dynOneEqEddy::correct(const tmp<volTensorField>& gradU)
{
    LESModel::correct(gradU);

    const volSymmTensorField D(symm(gradU));

    volScalarField KK(0.5*(filter_(magSqr(U())) - magSqr(filter_(U()))));
    KK.max(dimensionedScalar("small", KK.dimensions(), SMALL));

    const volScalarField P(2.0*nuSgs_*magSqr(D));

    tmp<fvScalarMatrix> kEqn
    (
       fvm::ddt(k_)
     + fvm::div(phi(), k_)
     - fvm::laplacian(DkEff(), k_)
    ==
       P
     - fvm::Sp(ce(D, KK)*sqrt(k_)/delta(), k_)
    );

    kEqn().relax();
    kEqn().solve();

    bound(k_, kMin_);

    updateSubGridScaleFields(D, KK);
}


bool dynOneEqEddy::read()
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
