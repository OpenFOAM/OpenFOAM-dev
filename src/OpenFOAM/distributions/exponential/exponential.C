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

#include "exponential.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace distributions
{
    defineTypeNameAndDebug(exponential, 0);
    addToRunTimeSelectionTable(distribution, exponential, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::distributions::exponential::phi
(
    const label q,
    const scalarField& x
) const
{
    return integerPow(x, q)*lambda_*exp(- lambda_*x);
}


Foam::tmp<Foam::scalarField> Foam::distributions::exponential::Phi
(
    const label q,
    const scalarField& x
) const
{
    if (q == 0)
    {
        return - exp(- lambda_*x);
    }
    else
    {
        return unintegrableForNonZeroQ::Phi(q, x);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distributions::exponential::exponential
(
    const unitConversion& units,
    const dictionary& dict,
    const label sampleQ,
    randomGenerator&& rndGen
)
:
    FieldDistribution<unintegrableForNonZeroQ, exponential>
    (
        typeName,
        units,
        dict,
        sampleQ,
        std::move(rndGen)
    ),
    min_(dict.lookupBackwardsCompatible<scalar>({"min", "minValue"}, units)),
    max_(dict.lookupBackwardsCompatible<scalar>({"max", "maxValue"}, units)),
    lambda_(dict.lookup<scalar>("lambda", unitless))
{
    validateBounds(dict);
    validatePositive(dict);
    mean();
}


Foam::distributions::exponential::exponential
(
    const exponential& d,
    const label sampleQ
)
:
    FieldDistribution<unintegrableForNonZeroQ, exponential>(d, sampleQ),
    min_(d.min_),
    max_(d.max_),
    lambda_(d.lambda_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::distributions::exponential::~exponential()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::distributions::exponential::sample() const
{
    if (q() == 0)
    {
        const scalar s = rndGen_.sample01<scalar>();
        const Pair<scalar>& Phi01 = this->Phi01();
        const scalar PhiS = (1 - s)*Phi01[0] + s*Phi01[1];
        return - 1/lambda_*log(- PhiS);
    }
    else
    {
        return unintegrableForNonZeroQ::sample();
    }
}


Foam::scalar Foam::distributions::exponential::min() const
{
    return min_;
}


Foam::scalar Foam::distributions::exponential::max() const
{
    return max_;
}


void Foam::distributions::exponential::write
(
    Ostream& os,
    const unitConversion& units
) const
{
    FieldDistribution<unintegrableForNonZeroQ, exponential>::write(os, units);

    writeEntry(os, "min", units, min_);
    writeEntry(os, "max", units, max_);
    writeEntry(os, "lambda", unitless, lambda_);
}


Foam::tmp<Foam::scalarField>
Foam::distributions::exponential::x(const label n) const
{
    tmp<scalarField> tx(distribution::x(n));
    tx.ref()[0] = Foam::max(tx.ref()[0], q() < 0 ? min_/2 : rootVSmall);
    return tx;
}


// ************************************************************************* //
