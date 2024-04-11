/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "normal.H"
#include "standardNormal.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace distributions
{
    defineTypeNameAndDebug(normal, 0);
    addToRunTimeSelectionTable(distribution, normal, dictionary);
}
}


using namespace Foam::constant::mathematical;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::distributions::normal::phi
(
    const label q,
    const scalarField& x
) const
{
    static const scalar sqrt2Pi = sqrt(2*pi);
    return integerPow(x, q)/(sigma_*sqrt2Pi)*exp(- sqr((x - mu_)/sigma_)/2);
}


Foam::tmp<Foam::scalarField> Foam::distributions::normal::Phi
(
    const label q,
    const scalarField& x
) const
{
    if (q == 0)
    {
        static const scalar sqrt2 = sqrt(scalar(2));
        return (1 + standardNormal::approxErf((x - mu_)/(sigma_*sqrt2)))/2;
    }
    else
    {
        return unintegrableForNonZeroQ::Phi(q, x);
    }
}


Foam::scalar Foam::distributions::normal::sampleForZeroQ(const scalar s) const
{
    static const scalar sqrt2 = sqrt(scalar(2));
    const Pair<scalar>& Phi01 = this->Phi01();
    const scalar PhiS = (1 - s)*Phi01[0] + s*Phi01[1];
    return standardNormal::approxErfInv(2*PhiS - 1)*sigma_*sqrt2 + mu_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distributions::normal::normal
(
    const dictionary& dict,
    randomGenerator& rndGen,
    const label sampleQ
)
:
    FieldDistribution<unintegrableForNonZeroQ, normal>
    (
        typeName,
        dict,
        rndGen,
        sampleQ
    ),
    min_(dict.lookupBackwardsCompatible<scalar>({"min", "minValue"})),
    max_(dict.lookupBackwardsCompatible<scalar>({"max", "maxValue"})),
    mu_(dict.lookupBackwardsCompatible<scalar>({"mu", "expectation"})),
    sigma_(dict.lookup<scalar>("sigma"))
{
    validateBounds(dict);
    if (q() != 0) validatePositive(dict);
    mean();
    report();
}


Foam::distributions::normal::normal
(
    randomGenerator& rndGen,
    const label Q,
    const label sampleQ,
    const label n,
    const scalar min,
    const scalar max,
    const scalar mu,
    const scalar sigma
)
:
    FieldDistribution<unintegrableForNonZeroQ, normal>(rndGen, Q, sampleQ, n),
    min_(min),
    max_(max),
    mu_(mu),
    sigma_(sigma)
{}


Foam::distributions::normal::normal(const normal& d, const label sampleQ)
:
    FieldDistribution<unintegrableForNonZeroQ, normal>(d, sampleQ),
    min_(d.min_),
    max_(d.max_),
    mu_(d.mu_),
    sigma_(d.sigma_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::distributions::normal::~normal()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::distributions::normal::sample() const
{
    if (q() == 0)
    {
        return sampleForZeroQ(rndGen_.sample01<scalar>());
    }
    else
    {
        return unintegrableForNonZeroQ::sample();
    }
}


Foam::scalar Foam::distributions::normal::min() const
{
    return min_;
}


Foam::scalar Foam::distributions::normal::max() const
{
    return max_;
}


Foam::scalar Foam::distributions::normal::mean() const
{
    if (q() == 0)
    {
        // !!! This isn't technically correct. The min/max clipping affects the
        // mean. To calculate it properly we have to integrate the first moment
        // of the PDF, and that can't be done analytically; we have to hand
        // over to the unintegrable handling. This probably isn't worth it for
        // what should be a small deviation. Also, it looks a bit odd in the
        // log for the distribution to be reporting a different mean to that
        // which was specified. So, just use the mean value of the un-clipped
        // distribution for now.
        return mu_;
    }
    else
    {
        return unintegrableForNonZeroQ::mean();
    }
}


Foam::tmp<Foam::scalarField>
Foam::distributions::normal::x(const label n) const
{
    tmp<scalarField> tx(distribution::x(n));
    tx.ref()[0] = Foam::max(tx.ref()[0], q() < 0 ? min_/2 : -vGreat);
    return tx;
}


// ************************************************************************* //
