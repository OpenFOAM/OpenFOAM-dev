/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2025 OpenFOAM Foundation
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

#include "fractal.H"
#include "populationBalanceModel.H"
#include "addToRunTimeSelectionTable.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalance
{
namespace shapeModels
{
    defineTypeNameAndDebug(fractal, 0);
    addToRunTimeSelectionTable(shapeModel, fractal, dictionary);
}
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::populationBalance::shapeModels::fractal::dColl(const label i) const
{
    using Foam::constant::mathematical::pi;

    tmp<volScalarField> tDColl
    (
        volScalarField::New
        (
            "dColl",
            popBal_.mesh(),
            dimensionedScalar(dimLength, Zero)
        )
    );

    volScalarField& dColl = tDColl.ref();

    dColl =
        6/kappas_[i]
       *pow(popBal_.v(i)*pow3(kappas_[i])/(36*pi*alphaC(i)), 1/Df(i));

    return tDColl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalance::shapeModels::fractal::fractal
(
    const dictionary& dict,
    const populationBalanceModel& popBal
)
:
    SecondaryPropertyModel<shapeModel>(popBal),
    Df_(Function1<scalar>::New("Df", {dimLength, dimless}, dict)),
    alphaC_(Function1<scalar>::New("alphaC", {dimLength, dimless}, dict)),
    kappas_(popBal.nGroups()),
    dColls_(popBal.nGroups()),
    Sus_(popBal.nGroups())
{
    forAll(popBal_.fs(), i)
    {
        kappas_.set
        (
            i,
            new volScalarField
            (
                populationBalanceModel::groupFieldIo
                (
                    "kappa",
                    i,
                    popBal_.phases()[i]
                ),
                populationBalanceModel::groupField
                (
                    "kappa",
                    i,
                    popBal_.phases()[i]
                )
            )
        );

        dColls_.set
        (
            i,
            new volScalarField
            (
                populationBalanceModel::groupFieldIo
                (
                    "dColl",
                    i,
                    popBal_.phases()[i]
                ),
                this->dColl(i)
            )
        );

        Sus_.set
        (
            i,
            new volScalarField::Internal
            (
                populationBalanceModel::groupFieldIo
                (
                    "kappa:Su",
                    i,
                    popBal_.phases()[i]
                ),
                popBal_.mesh(),
                dimensionedScalar(kappas_[i].dimensions()/dimTime, Zero)
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::populationBalance::shapeModels::fractal::~fractal()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::dimensionedScalar
Foam::populationBalance::shapeModels::fractal::Df(const label i) const
{
    return
        dimensionedScalar
        (
            "Df" + Foam::name(i),
            dimless,
            Df_->value(popBal_.dSph(i).value())
        );
}


Foam::dimensionedScalar
Foam::populationBalance::shapeModels::fractal::alphaC(const label i) const
{
    return
        dimensionedScalar
        (
            "alphaC" + Foam::name(i),
            dimless,
            alphaC_->value(popBal_.dSph(i).value())
        );
}


const Foam::volScalarField&
Foam::populationBalance::shapeModels::fractal::fld(const label i) const
{
    return kappas_[i];
}


Foam::volScalarField::Internal&
Foam::populationBalance::shapeModels::fractal::src(const label i)
{
    return Sus_[i];
}


void Foam::populationBalance::shapeModels::fractal::solve()
{
    forAll(popBal_.fs(), i)
    {
        const phaseModel& phase = popBal_.phases()[i];
        const volScalarField& alpha = phase;
        const volScalarField& rho = popBal_.phases()[i].rho();
        const volScalarField& fi = popBal_.fs()[i];

        const volScalarField alphaFi
        (
            IOobject::groupName
            (
                alpha.member() + fi.member().capitalise(),
                alpha.group()
            ),
            alpha*fi
        );

        const surfaceScalarField alphaFiPhi
        (
            IOobject::groupName
            (
                alpha.member() + fi.member().capitalise() + "Phi",
                alpha.group()
            ),
            max(fvc::interpolate(fi, "fi"), small)*phase.alphaPhi()
        );

        fvScalarMatrix kappaiEqn
        (
            fvm::ddt(alpha, fi, kappas_[i])
          + fvm::div(alphaFiPhi, kappas_[i])
         ==
            Sus_[i] + fvm::Sp(popBal_.Sp(i)*fi, kappas_[i])
          + popBal_.expansionSu(i, kappas_)
          + fvm::Sp(popBal_.expansionSp(i)*fi, kappas_[i])
          + popBal_.modelSourceSu(i, kappas_)
          + popBal_.fluid().fvModels().source(alphaFi, rho, kappas_[i])/rho
          - correction
            (
                fvm::Sp
                (
                    max(phase.residualAlpha() - alpha*fi, scalar(0))
                   /kappas_[i].mesh().time().deltaT(),
                    kappas_[i]
                )
            )
        );

        kappaiEqn.relax();

        popBal_.fluid().fvConstraints().constrain(kappaiEqn);

        kappaiEqn.solve();

        popBal_.fluid().fvConstraints().constrain(kappas_[i]);

        // Bound kappa so that the surface-area-volume ratio is greater than
        // that of spherical particles of this group, but less than that of the
        // particles represented by the first group
        kappas_[i] = min(max(kappas_[i], 6/popBal_.dSph(i)), 6/popBal_.dSph(0));

        kappas_[i].correctBoundaryConditions();
    }
}


void Foam::populationBalance::shapeModels::fractal::correct()
{
    forAll(popBal_.fs(), i)
    {
        // Update the collisional diameter
        dColls_[i] = dColl(i);
    }
}


Foam::tmp<Foam::volScalarField>
Foam::populationBalance::shapeModels::fractal::a(const label i) const
{
    return kappas_[i]*popBal_.v(i);
}


Foam::tmp<Foam::volScalarField>
Foam::populationBalance::shapeModels::fractal::d(const label i) const
{
    return dColls_[i];
}


// ************************************************************************* //
