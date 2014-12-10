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

#include "blackBodyEmission.H"
#include "physicoChemicalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::List<Foam::Tuple2<Foam::scalar, Foam::scalar> >
Foam::radiation::blackBodyEmission::emissivePowerTable
(
    IStringStream
    (
        "("
        "( 1000 0.00032)"
        "( 1100 0.00091)"
        "( 1200 0.00213)"
        "( 1300 0.00432)"
        "( 1400 0.00779)"
        "( 1500 0.01280)"
        "( 1600 0.01972)"
        "( 1700 0.02853)"
        "( 1800 0.03934)"
        "( 1900 0.05210)"
        "( 2000 0.06672)"
        "( 2100 0.08305)"
        "( 2200 0.10088)"
        "( 2300 0.12002)"
        "( 2400 0.14025)"
        "( 2500 0.16135)"
        "( 2600 0.18311)"
        "( 2700 0.20535)"
        "( 2800 0.22788)"
        "( 2900 0.25055)"
        "( 3000 0.27322)"
        "( 3100 0.29576)"
        "( 3200 0.31809)"
        "( 3300 0.34009)"
        "( 3400 0.36172)"
        "( 3500 0.38290)"
        "( 3600 0.40359)"
        "( 3700 0.42375)"
        "( 3800 0.44336)"
        "( 3900 0.46240)"
        "( 4000 0.48085)"
        "( 4100 0.49872)"
        "( 4200 0.51599)"
        "( 4300 0.53267)"
        "( 4400 0.54877)"
        "( 4500 0.56429)"
        "( 4600 0.57925)"
        "( 4700 0.59366)"
        "( 4800 0.60753)"
        "( 4900 0.62088)"
        "( 5000 0.63372)"
        "( 5100 0.64606)"
        "( 5200 0.65794)"
        "( 5300 0.66935)"
        "( 5400 0.68033)"
        "( 5500 0.69087)"
        "( 5600 0.70101)"
        "( 5700 0.71076)"
        "( 5800 0.72012)"
        "( 5900 0.72913)"
        "( 6000 0.73778)"
        "( 6100 0.74610)"
        "( 6200 0.75410)"
        "( 6300 0.76180)"
        "( 6400 0.76920)"
        "( 6500 0.77631)"
        "( 6600 0.78316)"
        "( 6700 0.78975)"
        "( 6800 0.79609)"
        "( 6900 0.80219)"
        "( 7000 0.80807)"
        "( 7100 0.81373)"
        "( 7200 0.81918)"
        "( 7300 0.82443)"
        "( 7400 0.82949)"
        "( 7500 0.83436)"
        "( 7600 0.83906)"
        "( 7700 0.84359)"
        "( 7800 0.84796)"
        "( 7900 0.85218)"
        "( 8000 0.85625)"
        "( 8100 0.86017)"
        "( 8200 0.86396)"
        "( 8300 0.86762)"
        "( 8400 0.87115)"
        "( 8500 0.87456)"
        "( 8600 0.87786)"
        "( 8700 0.88105)"
        "( 8800 0.88413)"
        "( 8900 0.88711)"
        "( 9000 0.88999)"
        "( 9100 0.89277)"
        "( 9200 0.89547)"
        "( 9300 0.89807)"
        "( 9400 0.90060)"
        "( 9500 0.90304)"
        "( 9600 0.90541)"
        "( 9700 0.90770)"
        "( 9800 0.90992)"
        "( 9900 0.91207)"
        "(10000 0.91415)"
        "(12000 0.94505)"
        "(15000 0.96893)"
        "(20000 0.98555)"
        "(30000 0.99529)"
        "(40000 0.99792)"
        "(50000 0.99890)"
        ")"
    )()
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::blackBodyEmission::blackBodyEmission
(
    const label nLambda,
    const volScalarField& T
)
:
    table_
    (
        emissivePowerTable,
        interpolationTable<scalar>::CLAMP,
        "blackBodyEmissivePower"
    ),
    C1_("C1", dimensionSet(1, 4, 3, 0, 0, 0, 0), 3.7419e-16),
    C2_("C2", dimensionSet(0, 1, 0, 1, 0, 0, 0), 14.388e-6),
    bLambda_(nLambda),
    T_(T)
{
    forAll(bLambda_, lambdaI)
    {
        bLambda_.set
        (
            lambdaI,
            new volScalarField
            (
                IOobject
                (
                    "bLambda_" + Foam::name(lambdaI) ,
                    T.mesh().time().timeName(),
                    T.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                physicoChemical::sigma*pow4(T)
            )
        );

    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::blackBodyEmission::~blackBodyEmission()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::radiation::blackBodyEmission::fLambdaT
(
    const scalar lambdaT
) const
{
    return  table_(lambdaT*1.0e6);
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::blackBodyEmission::EbDeltaLambdaT
(
    const volScalarField& T,
    const Vector2D<scalar>& band
) const
{
    tmp<volScalarField> Eb
    (
        new volScalarField
        (
            IOobject
            (
                "Eb",
                T.mesh().time().timeName(),
                T.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            physicoChemical::sigma*pow4(T)
        )
    );


    if (band == Vector2D<scalar>::one)
    {
        return Eb;
    }
    else
    {
        forAll(T, i)
        {
            scalar T1 = fLambdaT(band[1]*T[i]);
            scalar T2 = fLambdaT(band[0]*T[i]);
            dimensionedScalar fLambdaDelta
            (
                "fLambdaDelta",
                dimless,
                T1 - T2
            );
            Eb()[i] = Eb()[i]*fLambdaDelta.value();
        }
        return Eb;
    }
}


void Foam::radiation::blackBodyEmission::correct
(
    const label lambdaI,
    const Vector2D<scalar>& band
)
{
    bLambda_[lambdaI] = EbDeltaLambdaT(T_, band);
}


// ************************************************************************* //
