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

#include "RosinRammler.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace distributions
{
    defineTypeNameAndDebug(RosinRammler, 0);
    addToRunTimeSelectionTable(distribution, RosinRammler, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::distributions::RosinRammler::phi
(
    const label q,
    const scalarField& x
) const
{
    const scalarField xByDPowNm1(pow(x/d_, n_ - 1));

    return integerPow(x, q)*(n_/d_)*xByDPowNm1*exp(- xByDPowNm1*x/d_);
}


Foam::tmp<Foam::scalarField> Foam::distributions::RosinRammler::Phi
(
    const label q,
    const scalarField& x
) const
{
    if (q == 0)
    {
        return - exp(- pow(x/d_, n_));
    }
    else
    {
        return unintegrableForNonZeroQ::Phi(q, x);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distributions::RosinRammler::RosinRammler
(
    const unitConversion& units,
    const dictionary& dict,
    const label sampleQ,
    randomGenerator&& rndGen
)
:
    FieldDistribution<unintegrableForNonZeroQ, RosinRammler>
    (
        typeName,
        units,
        dict,
        sampleQ,
        std::move(rndGen)
    ),
    min_(dict.lookupBackwardsCompatible<scalar>({"min", "minValue"}, units)),
    max_(dict.lookupBackwardsCompatible<scalar>({"max", "maxValue"}, units)),
    d_(dict.lookup<scalar>("d", units)),
    n_(dict.lookup<scalar>("n", unitless))
{
    validateBounds(dict);
    validatePositive(dict);
    mean();
}


Foam::distributions::RosinRammler::RosinRammler
(
    const RosinRammler& d,
    const label sampleQ
)
:
    FieldDistribution<unintegrableForNonZeroQ, RosinRammler>(d, sampleQ),
    min_(d.min_),
    max_(d.max_),
    d_(d.d_),
    n_(d.n_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::distributions::RosinRammler::~RosinRammler()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::distributions::RosinRammler::sample() const
{
    if (q() == 0)
    {
        const scalar s = rndGen_.sample01<scalar>();
        const Pair<scalar>& Phi01 = this->Phi01();
        const scalar PhiS = (1 - s)*Phi01[0] + s*Phi01[1];
        return d_*pow(- log(- PhiS), 1/n_);
    }
    else
    {
        return unintegrableForNonZeroQ::sample();
    }
}


Foam::scalar Foam::distributions::RosinRammler::min() const
{
    return min_;
}


Foam::scalar Foam::distributions::RosinRammler::max() const
{
    return max_;
}


void Foam::distributions::RosinRammler::write
(
    Ostream& os,
    const unitConversion& units
) const
{
    FieldDistribution<unintegrableForNonZeroQ, RosinRammler>::write(os, units);

    writeEntry(os, "min", units, min_);
    writeEntry(os, "max", units, max_);
    writeEntry(os, "d", units, d_);
    writeEntry(os, "n", unitless, n_);
}


Foam::tmp<Foam::scalarField>
Foam::distributions::RosinRammler::x(const label n) const
{
    tmp<scalarField> tx(distribution::x(n));
    tx.ref()[0] = Foam::max(tx.ref()[0], q() < 0 ? min_/2 : rootVSmall);
    return tx;
}


// ************************************************************************* //
