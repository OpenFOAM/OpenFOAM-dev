/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2026 OpenFOAM Foundation
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
    const dictionary& coeffDict,
    const uRhoMulticomponentThermo& ct
)
:
    laminarFlameSpeed(dict, ct),
    pPoints_(coeffDict.lookup("pPoints")),
    PhiPoints_(coeffDict.lookup("PhiPoints")),
    alpha_(coeffDict.lookup("alpha")),
    beta_(coeffDict.lookup("beta")),
    TRef_(coeffDict.lookup<scalar>("TRef"))
{
    checkPointsMonotonicity(coeffDict, "Phi", PhiPoints_);
    checkPointsMonotonicity(coeffDict, "pressure", pPoints_);
    checkCoefficientArrayShape(coeffDict, "alpha", alpha_);
    checkCoefficientArrayShape(coeffDict, "beta", beta_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::laminarFlameSpeedModels::RaviPetersen::~RaviPetersen()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::laminarFlameSpeedModels::RaviPetersen::checkPointsMonotonicity
(
    const dictionary& coeffDict,
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
                coeffDict
            )   << "Data points for the " << name
                << " do not increase monotonically" << endl
                << exit(FatalIOError);
        }
    }
}


void Foam::laminarFlameSpeedModels::RaviPetersen::checkCoefficientArrayShape
(
    const dictionary& coeffDict,
    const word& name,
    const List<List<List<scalar>>>& x
) const
{
    bool ok = true;

    ok &= x.size() == PhiPoints_.size() - 1;

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
            coeffDict
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
    const label PhiIndex,
    const label pIndex,
    const scalar Phi,
    const scalar Tu
) const
{
    return pow
    (
        Tu/TRef_,
        polynomial(beta_[PhiIndex][pIndex],Phi)
    );
}


inline Foam::scalar
Foam::laminarFlameSpeedModels::RaviPetersen::correlationInRange
(
    const label PhiIndex,
    const label pIndex,
    const scalar Phi,
    const scalar Tu
) const
{
    // standard correlation
    return
        polynomial(alpha_[PhiIndex][pIndex],Phi)
       *THatPowB(PhiIndex, pIndex, Phi, Tu);
}


inline Foam::scalar
Foam::laminarFlameSpeedModels::RaviPetersen::correlationOutOfRange
(
    const label PhiIndex,
    const label pIndex,
    const scalar Phi,
    const scalar PhiLim,
    const scalar Tu
) const
{
    scalar A = polynomial(alpha_[PhiIndex][pIndex], PhiLim);
    scalar dA = dPolynomial(alpha_[PhiIndex][pIndex], PhiLim);
    scalar dB = dPolynomial(beta_[PhiIndex][pIndex], PhiLim);
    scalar TB = THatPowB(PhiIndex, pIndex, PhiLim, Tu);

    // linear extrapolation from the bounds of the correlation
    return max(TB*(A + (dA + A*log(Tu/TRef_)*dB)*(Phi - PhiLim)), 0.0);
}


inline Foam::scalar Foam::laminarFlameSpeedModels::RaviPetersen::Su0pTphi
(
    const scalar p,
    const scalar Tu,
    const scalar Phi
) const
{
    scalar Su = 0, s;

    label PhiIndex, pIndex;
    scalar PhiXi, pXi;
    scalar PhiLim, pLim;
    bool PhiInRange;

    PhiInRange = interval(PhiPoints_, Phi, PhiIndex, PhiXi, PhiLim);

    interval(pPoints_, p, pIndex, pXi, pLim);

    for (label pI = 0; pI < 2; pI ++)
    {
        if (PhiInRange)
        {
            s = correlationInRange(PhiIndex, pIndex + pI, Phi, Tu);
        }
        else
        {
            s = correlationOutOfRange(PhiIndex, pIndex + pI, Phi, PhiLim, Tu);
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
    if (uThermo_.containsSpecie("egr"))
    {
        FatalErrorInFunction
            << "The " << type() << " model does not support EGR"
            << exit(FatalError);
    }

    const volScalarField& p = uThermo_.p();
    const volScalarField& Tu = uThermo_.T();
    const volScalarField Phi(uThermo_.Phi());

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
        Su0[celli] = Su0pTphi(p[celli], Tu[celli], Phi[celli]);
    }

    volScalarField::Boundary& Su0Bf = Su0.boundaryFieldRef();

    forAll(Su0Bf, patchi)
    {
        forAll(Su0Bf[patchi], facei)
        {
            Su0Bf[patchi][facei] =
                Su0pTphi
                (
                    p.boundaryField()[patchi][facei],
                    Tu.boundaryField()[patchi][facei],
                    Phi.boundaryField()[patchi][facei]
                );
        }
    }

    return tSu0;
}


// ************************************************************************* //
