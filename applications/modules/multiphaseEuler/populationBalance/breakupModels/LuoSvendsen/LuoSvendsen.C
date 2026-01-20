/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2026 OpenFOAM Foundation
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

#include "LuoSvendsen.H"
#include "phaseCompressibleMomentumTransportModel.H"
#include "linearInterpolationWeights.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalance
{
namespace breakupModels
{
    defineTypeNameAndDebug(LuoSvendsen, 0);
    addToRunTimeSelectionTable(breakupModel, LuoSvendsen, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalance::breakupModels::LuoSvendsen::LuoSvendsen
(
    const populationBalanceModel& popBal,
    const dictionary& dict
)
:
    binary(popBal, dict),
    gammaUpperReg2by11_(),
    gammaUpperReg5by11_(),
    gammaUpperReg8by11_(),
    C4_("C4", dimless, dict, 0.923),
    beta_("beta", dimless, dict, 2.05),
    minEddyRatio_("minEddyRatio", dimless, dict, 11.4),
    kolmogorovLengthScale_
    (
        IOobject
        (
            "kolmogorovLengthScale",
            popBal_.time().name(),
            popBal_.mesh()
        ),
        popBal_.mesh(),
        dimensionedScalar
        (
            "kolmogorovLengthScale",
            dimLength,
            Zero
        )
    )
{
    List<Tuple2<scalar, scalar>> gammaUpperReg2by11Table;
    List<Tuple2<scalar, scalar>> gammaUpperReg5by11Table;
    List<Tuple2<scalar, scalar>> gammaUpperReg8by11Table;

    gammaUpperReg2by11Table.append(Tuple2<scalar, scalar>(0.0, 1.0));
    gammaUpperReg5by11Table.append(Tuple2<scalar, scalar>(0.0, 1.0));
    gammaUpperReg8by11Table.append(Tuple2<scalar, scalar>(0.0, 1.0));

    for (scalar z = 1e-2; z <= 10.0; z = z + 1e-2)
    {
        Tuple2<scalar, scalar> gamma2by11
            (
                z,
                incGammaRatio_Q(2.0/11.0, z)
            );

        Tuple2<scalar, scalar> gamma5by11
            (
                z,
                incGammaRatio_Q(5.0/11.0, z)
            );

        Tuple2<scalar, scalar> gamma8by11
            (
                z,
                incGammaRatio_Q(8.0/11.0, z)
            );

        gammaUpperReg2by11Table.append(gamma2by11);
        gammaUpperReg5by11Table.append(gamma5by11);
        gammaUpperReg8by11Table.append(gamma8by11);
    }

    gammaUpperReg2by11_ =
        new Function1s::Table<scalar>
        (
            "gamma2by11",
            Function1s::tableBase::boundsHandling::clamp,
            linearInterpolationWeights::typeName,
            autoPtr<TableReader<scalar, scalar>>(nullptr),
            gammaUpperReg2by11Table
        );

    gammaUpperReg5by11_ =
        new Function1s::Table<scalar>
        (
            "gamma5by11",
            Function1s::tableBase::boundsHandling::clamp,
            linearInterpolationWeights::typeName,
            autoPtr<TableReader<scalar, scalar>>(nullptr),
            gammaUpperReg5by11Table
        );

    gammaUpperReg8by11_ =
        new Function1s::Table<scalar>
        (
            "gamma8by11",
            Function1s::tableBase::boundsHandling::clamp,
            linearInterpolationWeights::typeName,
            autoPtr<TableReader<scalar, scalar>>(nullptr),
            gammaUpperReg8by11Table
        );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::populationBalance::breakupModels::LuoSvendsen::precompute()
{
    kolmogorovLengthScale_ =
        pow025
        (
            pow3
            (
                popBal_.continuousPhase().fluidThermo().nu()
            )
           /popBal_.continuousTurbulence().epsilon()
        );
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::populationBalance::breakupModels::LuoSvendsen::rate
(
    const label i,
    const label j
) const
{
    const dimensionedScalar& dSphj = popBal_.dSph(j);
    const dimensionedScalar& vi = popBal_.v(i);
    const dimensionedScalar& vj = popBal_.v(j);

    const volScalarField::Internal& rhoc = popBal_.continuousPhase().rho();

    tmp<volScalarField> tsigma(popBal_.sigmaWithContinuousPhase(i));
    const volScalarField::Internal& sigma = tsigma();

    tmp<volScalarField> tepsilonc(popBal_.continuousTurbulence().epsilon());
    const volScalarField::Internal& epsilonc = tepsilonc();

    const dimensionedScalar cf
    (
        pow(vi/vj, 2.0/3.0) + pow((1 - vi/vj), 2.0/3.0) - 1
    );

    const volScalarField::Internal b
    (
        12*cf*sigma
       /(beta_*rhoc*pow(dSphj, 5.0/3.0)*pow(epsilonc, 2.0/3.0))
    );

    const volScalarField::Internal xiMin
    (
        minEddyRatio_*kolmogorovLengthScale_/dSphj
    );

    const volScalarField::Internal tMin(b/pow(xiMin, 11.0/3.0));

    volScalarField::Internal integral(3/(11*pow(b, 8.0/11.0)));
    forAll(integral, celli)
    {
        integral[celli] *=
            2
           *pow(b[celli], 3.0/11.0)
           *tgamma(5.0/11.0)
           *(
                gammaUpperReg5by11_->value(b[celli])
              - gammaUpperReg5by11_->value(tMin[celli])
            );
    }

    return
        C4_
       *(1 - popBal_.alphas()())
       /vj
       *cbrt(epsilonc/sqr(dSphj))
       *integral;
}


// ************************************************************************* //
