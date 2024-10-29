/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "standardNormal.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace distributions
{
    defineTypeNameAndDebug(standardNormal, 0);
}
}


const Foam::scalar Foam::distributions::standardNormal::a_ = 0.147;


using namespace Foam::constant::mathematical;


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::distributions::standardNormal::approxErf
(
    const scalarField& x
)
{
    return sign(x)*sqrt(1 - exp(-sqr(x)*(4/pi + a_*sqr(x))/(1 + a_*sqr(x))));
}


Foam::scalar Foam::distributions::standardNormal::approxErfInv(const scalar y)
{
    const scalar l = log(Foam::max(1 - y*y, small/2)), b = 2/(pi*a_) + l/2;
    return sign(y)*sqrt(-b + sqrt(b*b - l/a_));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distributions::standardNormal::standardNormal(randomGenerator&& rndGen)
:
    FieldDistribution<distribution, standardNormal>(0, 0, std::move(rndGen))
{}


Foam::distributions::standardNormal::standardNormal
(
    const randomGenerator::seed& s,
    const bool global
)
:
    standardNormal(randomGenerator(s, global))
{}


Foam::distributions::standardNormal::standardNormal(const standardNormal& d)
:
    FieldDistribution<distribution, standardNormal>(d, 0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::distributions::standardNormal::~standardNormal()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::distributions::standardNormal::sample() const
{
    static const scalar sqrt2 = sqrt(scalar(2));
    const scalar s = rndGen_.sample01<scalar>();
    return approxErfInv(2*s - 1)*sqrt2;
}


Foam::scalar Foam::distributions::standardNormal::min() const
{
    return -vGreat;
}


Foam::scalar Foam::distributions::standardNormal::max() const
{
    return vGreat;
}


Foam::scalar Foam::distributions::standardNormal::mean() const
{
    return 0;
}


Foam::tmp<Foam::scalarField>
Foam::distributions::standardNormal::integralPDFxPow
(
    const scalarField& x,
    const label e,
    const bool
) const
{
    if (e != 0) NotImplemented;

    static const scalar sqrt2 = sqrt(scalar(2));
    return (1 + approxErf(x/sqrt2))/2;
}


Foam::tmp<Foam::scalarField>
Foam::distributions::standardNormal::plotX(const label n) const
{
    const scalar x0 = approxErfInv(1 - rootSmall);
    const scalar x1 = approxErfInv(rootSmall - 1);

    tmp<scalarField> tResult(new scalarField(n));
    scalarField& result = tResult.ref();

    forAll(result, i)
    {
        const scalar f = scalar(i)/(n - 1);
        result[i] = (1 - f)*x0 + f*x1;
    }

    return tResult;
}


Foam::tmp<Foam::scalarField> Foam::distributions::standardNormal::plotPDF
(
    const scalarField& x
) const
{
    static const scalar sqrt2Pi = sqrt(2*pi);
    return exp(- sqr(x)/2)/sqrt2Pi;
}


// ************************************************************************* //
