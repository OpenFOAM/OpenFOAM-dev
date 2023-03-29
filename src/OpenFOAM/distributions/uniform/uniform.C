/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "uniform.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace distributions
{
    defineTypeNameAndDebug(uniform, 0);
    addToRunTimeSelectionTable(distribution, uniform, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::distributions::uniform::Phi(const scalar x) const
{
    if (q() == -1)
    {
        return log(x);
    }
    else
    {
        return integerPow(x, 1 + q())/(1 + q());
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distributions::uniform::uniform
(
    const dictionary& dict,
    Random& rndGen,
    const label sampleQ
)
:
    FieldDistribution<distribution, uniform>(typeName, dict, rndGen, sampleQ),
    min_(dict.lookupBackwardsCompatible<scalar>({"min", "minValue"})),
    max_(dict.lookupBackwardsCompatible<scalar>({"max", "maxValue"})),
    Phi0_(Phi(min_)),
    Phi1_(Phi(max_))
{
    validateBounds(dict);
    if (q() != 0) validatePositive(dict);
    report();
}


Foam::distributions::uniform::uniform(const uniform& d, const label sampleQ)
:
    FieldDistribution<distribution, uniform>(d, sampleQ),
    min_(d.min_),
    max_(d.max_),
    Phi0_(Phi(min_)),
    Phi1_(Phi(max_))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::distributions::uniform::~uniform()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::distributions::uniform::sample() const
{
    const scalar s = rndGen_.sample01<scalar>();

    if (q() == -1)
    {
        return min_*pow(max_/min_, s);
    }
    else
    {
        const scalar PhiS = (1 - s)*Phi0_ + s*Phi1_;

        return integerRoot((1 + q())*PhiS, 1 + q());
    }
}


Foam::scalar Foam::distributions::uniform::min() const
{
    return min_;
}


Foam::scalar Foam::distributions::uniform::max() const
{
    return max_;
}


Foam::scalar Foam::distributions::uniform::mean() const
{
    if (q() == -2)
    {
        return Foam::log(max_/min_)/(Phi1_ - Phi0_);
    }
    else
    {
        const scalar Mu0 = integerPow(min_, 2 + q())/(2 + q());
        const scalar Mu1 = integerPow(max_, 2 + q())/(2 + q());

        return (Mu1 - Mu0)/(Phi1_ - Phi0_);
    }
}


Foam::tmp<Foam::scalarField>
Foam::distributions::uniform::PDF(const scalarField& x) const
{
    if (q() == -1)
    {
        return clipPDF(x, 1/x/Foam::log(max_/min_));
    }
    else
    {
        return clipPDF(x, integerPow(x, q())/(Phi1_ - Phi0_));
    }
}


// ************************************************************************* //
