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

#include "RaviPetersen.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarFlameSpeedModels
{
    defineTypeNameAndDebug(RaviPetersen, 0);

    addToRunTimeSelectionTable
    (
        laminarFlameSpeed,
        RaviPetersen,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::laminarFlameSpeedModels::RaviPetersen::RaviPetersen
(
    const dictionary& dict,
    const psiuMulticomponentThermo& ct
)
:
    laminarFlameSpeed(dict, ct),
    coeffsDict_(dict.optionalSubDict(typeName + "Coeffs").subDict(fuel_)),
    pPoints_(coeffsDict_.lookup("pPoints")),
    EqRPoints_(coeffsDict_.lookup("EqRPoints")),
    alpha_(coeffsDict_.lookup("alpha")),
    beta_(coeffsDict_.lookup("beta")),
    TRef_(coeffsDict_.lookup<scalar>("TRef"))
{
    checkPointsMonotonicity("equivalenceRatio", EqRPoints_);
    checkPointsMonotonicity("pressure", pPoints_);
    checkCoefficientArrayShape("alpha", alpha_);
    checkCoefficientArrayShape("beta", beta_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::laminarFlameSpeedModels::RaviPetersen::~RaviPetersen()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::laminarFlameSpeedModels::RaviPetersen::checkPointsMonotonicity
(
    const word& name,
    const List<scalar>& x
) const
{
    for (label i = 1; i < x.size(); i ++)
    {
        if (x[i] <= x[i-1])
        {
            FatalIOErrorInFunction
            (
                coeffsDict_
            )   << "Data points for the " << name
                << " do not increase monotonically" << endl
                << exit(FatalIOError);
        }
    }
}


void Foam::laminarFlameSpeedModels::RaviPetersen::checkCoefficientArrayShape
(
    const word& name,
    const List<List<List<scalar>>>& x
) const
{
    bool ok = true;

    ok &= x.size() == EqRPoints_.size() - 1;

    forAll(x, i)
    {
        ok &= x[i].size() == pPoints_.size();

        forAll(x[i], j)
        {
            ok &= x[i][j].size() == x[i][0].size();
        }
    }

    if (!ok)
    {
        FatalIOErrorInFunction
        (
            coeffsDict_
        )   << "Inconsistent size of " << name << " coefficients array" << endl
            << exit(FatalIOError);
    }
}


inline bool Foam::laminarFlameSpeedModels::RaviPetersen::interval
(
    const List<scalar>& xPoints,
    const scalar x,
    label& xIndex,
    scalar& xXi,
    scalar& xLim
) const
{
    if (x < xPoints.first())
    {
        xIndex = 0;
        xXi = 0.0;
        xLim = xPoints.first();
        return false;
    }

    else if (x > xPoints.last())
    {
        xIndex = xPoints.size() - 2;
        xXi = 1.0;
        xLim = xPoints.last();
        return false;
    }

    for (xIndex = 0; x > xPoints[xIndex+1]; xIndex ++)
    {
        // increment xIndex until xPoints[xIndex] < x < xPoints[xIndex+1]
    }

    xXi = (x - xPoints[xIndex])/(xPoints[xIndex+1] - xPoints[xIndex]);
    xLim = x;

    return true;
}


inline Foam::scalar Foam::laminarFlameSpeedModels::RaviPetersen::polynomial
(
    const List<scalar>& coeffs,
    const scalar x
) const
{
    scalar xPow = 1.0;
    scalar y = 0.0;
    forAll(coeffs, i)
    {
        y += coeffs[i]*xPow;
        xPow *= x;
    }
    return y;
}


inline Foam::scalar Foam::laminarFlameSpeedModels::RaviPetersen::dPolynomial
(
    const List<scalar>& coeffs,
    const scalar x
) const
{
    scalar xPow = 1.0;
    scalar y = 0.0;
    for (label i = 1; i < coeffs.size(); i++)
    {
        y += i*coeffs[i]*xPow;
        xPow *= x;
    }
    return y;
}


inline Foam::scalar Foam::laminarFlameSpeedModels::RaviPetersen::THatPowB
(
    const label EqRIndex,
    const label pIndex,
    const scalar EqR,
    const scalar Tu
) const
{
    return pow
    (
        Tu/TRef_,
        polynomial(beta_[EqRIndex][pIndex],EqR)
    );
}


inline Foam::scalar
Foam::laminarFlameSpeedModels::RaviPetersen::correlationInRange
(
    const label EqRIndex,
    const label pIndex,
    const scalar EqR,
    const scalar Tu
) const
{
    // standard correlation
    return
        polynomial(alpha_[EqRIndex][pIndex],EqR)
       *THatPowB(EqRIndex, pIndex, EqR, Tu);
}


inline Foam::scalar
Foam::laminarFlameSpeedModels::RaviPetersen::correlationOutOfRange
(
    const label EqRIndex,
    const label pIndex,
    const scalar EqR,
    const scalar EqRLim,
    const scalar Tu
) const
{
    scalar A = polynomial(alpha_[EqRIndex][pIndex], EqRLim);
    scalar dA = dPolynomial(alpha_[EqRIndex][pIndex], EqRLim);
    scalar dB = dPolynomial(beta_[EqRIndex][pIndex], EqRLim);
    scalar TB = THatPowB(EqRIndex, pIndex, EqRLim, Tu);

    // linear extrapolation from the bounds of the correlation
    return max(TB*(A + (dA + A*log(Tu/TRef_)*dB)*(EqR - EqRLim)), 0.0);
}


inline Foam::scalar Foam::laminarFlameSpeedModels::RaviPetersen::speed
(
    const scalar EqR,
    const scalar p,
    const scalar Tu
) const
{
    scalar Su = 0, s;

    label EqRIndex, pIndex;
    scalar EqRXi, pXi;
    scalar EqRLim, pLim;
    bool EqRInRange;

    EqRInRange = interval(EqRPoints_, EqR, EqRIndex, EqRXi, EqRLim);

    interval(pPoints_, p, pIndex, pXi, pLim);

    for (label pI = 0; pI < 2; pI ++)
    {
        if (EqRInRange)
        {
            s = correlationInRange(EqRIndex, pIndex + pI, EqR, Tu);
        }
        else
        {
            s = correlationOutOfRange(EqRIndex, pIndex + pI, EqR, EqRLim, Tu);
        }

        Su += (1 - pXi)*s;
        pXi = 1 - pXi;
    }

    return Su;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::laminarFlameSpeedModels::RaviPetersen::operator()() const
{
    const volScalarField& p = psiuMulticomponentThermo_.p();
    const volScalarField& Tu = psiuMulticomponentThermo_.Tu();

    volScalarField EqR
    (
        IOobject
        (
            "EqR",
            p.time().name(),
            p.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        p.mesh(),
        dimensionedScalar(dimless, 0)
    );

    if (psiuMulticomponentThermo_.composition().contains("ft"))
    {
        const volScalarField& ft =
            psiuMulticomponentThermo_.composition().Y("ft");

        EqR =
            dimensionedScalar
            (
                "stoichiometricAirFuelMassRatio",
                dimless,
                psiuMulticomponentThermo_.properties()
            )*ft/max(1 - ft, small);
    }
    else
    {
        EqR = equivalenceRatio_;
    }

    tmp<volScalarField> tSu0
    (
        volScalarField::New
        (
            "Su0",
            p.mesh(),
            dimensionedScalar(dimVelocity, 0)
        )
    );

    volScalarField& Su0 = tSu0.ref();

    forAll(Su0, celli)
    {
        Su0[celli] = speed(EqR[celli], p[celli], Tu[celli]);
    }

    return tSu0;
}


// ************************************************************************* //
