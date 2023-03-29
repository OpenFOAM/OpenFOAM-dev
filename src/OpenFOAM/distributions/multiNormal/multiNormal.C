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

#include "multiNormal.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace distributions
{
    defineTypeNameAndDebug(multiNormal, 0);
    addToRunTimeSelectionTable(distribution, multiNormal, dictionary);
}
}


// * * * * * * * * * * * Private Static Member Functions * * * * * * * * * * //

Foam::scalarList Foam::distributions::multiNormal::readCumulativeStrengths
(
    const dictionary& dict
)
{
    const scalarList s(dict.lookup<scalarList>("strength"));

    const scalarList sHat(s/sum(s));

    scalarList cumSHat(s.size() + 1, scalar(0));
    forAll(sHat, i)
    {
        cumSHat[i + 1] = cumSHat[i] + sHat[i];
    }

    return cumSHat;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::distributions::multiNormal::phi
(
    const label q,
    const scalarField& x
) const
{
    scalarField phiQ0(x.size(), 0);

    forAll(distributions_, i)
    {
        phiQ0 +=
            (cumulativeStrengths_[i + 1] - cumulativeStrengths_[i])
           *distributions_[i].phi(0, x);
    }

    return integerPow(x, q)*phiQ0;
}


Foam::tmp<Foam::scalarField> Foam::distributions::multiNormal::Phi
(
    const label q,
    const scalarField& x
) const
{
    if (q == 0)
    {
        tmp<scalarField> tPhiQ0(new scalarField(x.size(), 0));
        scalarField& PhiQ0 = tPhiQ0.ref();

        forAll(distributions_, i)
        {
            PhiQ0 +=
                (cumulativeStrengths_[i + 1] - cumulativeStrengths_[i])
               *distributions_[i].Phi(0, x);
        }

        return tPhiQ0;
    }
    else
    {
        return unintegrableForNonZeroQ::Phi(q, x);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distributions::multiNormal::multiNormal
(
    const dictionary& dict,
    Random& rndGen,
    const label sampleQ
)
:
    FieldDistribution<unintegrableForNonZeroQ, multiNormal>
    (
        typeName,
        dict,
        rndGen,
        sampleQ
    ),
    cumulativeStrengths_(readCumulativeStrengths(dict))
{
    const scalar min
    (
        dict.lookupBackwardsCompatible<scalar>({"min", "minValue"})
    );
    const scalar max
    (
        dict.lookupBackwardsCompatible<scalar>({"max", "maxValue"})
    );
    const scalarList mu
    (
        dict.lookupBackwardsCompatible<scalarList>({"mu", "expectation"})
    );
    const scalarList sigma
    (
        dict.lookup<scalarList>("sigma")
    );

    if
    (
        mu.size() != cumulativeStrengths_.size() - 1
     || sigma.size() != cumulativeStrengths_.size() - 1
    )
    {
        FatalIOErrorInFunction(dict)
            << type() << ": Differing numbers of means, standard deviations "
            << "and strengths were given:" << nl << "    mu = " << mu
            << ", sigma = " << sigma << ", strength = "
            << dict.lookup<scalarList>("strength") << abort(FatalIOError);
    }

    distributions_.resize(mu.size());
    forAll(mu, i)
    {
        distributions_.set
        (
            i,
            new normal
            (
                rndGen_,
                0,
                0,
                -1,
                min,
                max,
                mu[i],
                sigma[i]
            )
        );
    }

    validateBounds(dict);
    if (q() != 0) validatePositive(dict);
    mean();
    report();
}


Foam::distributions::multiNormal::multiNormal
(
    const multiNormal& d,
    const label sampleQ
)
:
    FieldDistribution<unintegrableForNonZeroQ, multiNormal>(d, sampleQ),
    cumulativeStrengths_(d.cumulativeStrengths_),
    distributions_(d.distributions_.size())
{
    forAll(d.distributions_, i)
    {
        distributions_.set
        (
            i,
            new normal
            (
                rndGen_,
                0,
                0,
                -1,
                d.distributions_[i].min_,
                d.distributions_[i].max_,
                d.distributions_[i].mu_,
                d.distributions_[i].sigma_
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::distributions::multiNormal::~multiNormal()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::distributions::multiNormal::sample() const
{
    if (q() == 0)
    {
        const scalar S = rndGen_.sample01<scalar>();

        const label n = cumulativeStrengths_.size() - 1;
        label i = 0;
        for (; i < n && cumulativeStrengths_[i + 1] < S; ++ i);

        const scalar s =
            (S - cumulativeStrengths_[i])
           /(cumulativeStrengths_[i + 1] - cumulativeStrengths_[i]);

        return distributions_[i].sampleForZeroQ(s);
    }
    else
    {
        return unintegrableForNonZeroQ::sample();
    }
}


Foam::scalar Foam::distributions::multiNormal::min() const
{
    return distributions_[0].min();
}


Foam::scalar Foam::distributions::multiNormal::max() const
{
    return distributions_[0].max();
}


Foam::tmp<Foam::scalarField>
Foam::distributions::multiNormal::x(const label n) const
{
    return distributions_[0].x(n);
}


// ************************************************************************* //
