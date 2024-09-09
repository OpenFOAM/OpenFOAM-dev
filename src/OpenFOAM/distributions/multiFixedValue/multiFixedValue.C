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

#include "multiFixedValue.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace distributions
{
    defineTypeNameAndDebug(multiFixedValue, 0);
    addToRunTimeSelectionTable(distribution, multiFixedValue, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distributions::multiFixedValue::multiFixedValue
(
    const unitConversion& defaultUnits,
    const dictionary& dict,
    const label sampleQ,
    randomGenerator&& rndGen
)
:
    FieldDistribution<distribution, multiFixedValue>
    (
        typeName,
        defaultUnits,
        dict,
        sampleQ,
        std::move(rndGen)
    )
{
    List<Tuple2<scalar, scalar>> values(dict.lookup("values"));

    // Sort
    Foam::sort
    (
        values,
        [](const Tuple2<scalar, scalar>& a, const Tuple2<scalar, scalar>& b)
        {
            return a.first() < b.first();
        }
    );

    // Checks
    for (label i = 1; i < values.size(); ++ i)
    {
        if (values[i].second() < 0)
        {
            FatalIOErrorInFunction(dict)
                << typeName << ": The probabilities are not all positive "
                << abort(FatalIOError);
        }
    }

    // Optionally read units
    unitConversion units(defaultUnits);
    units.readIfPresent("units", dict);

    // Copy the coordinates
    x_.resize(values.size());
    forAll(values, i)
    {
        x_[i] = units.toStandard(values[i].first());
    }

    // Copy the probabilities. Scale if q != 0.
    P_.resize(values.size());
    forAll(values, i)
    {
        P_[i] = integerPow(x_[i], q())*values[i].second();
    }

    // Compute the cumulative sum of the probabilities
    sumP_.resize(values.size() + 1);
    sumP_[0] = 0;
    forAll(values, i)
    {
        sumP_[i + 1] = sumP_[i] + P_[i];
    }

    // Normalise
    P_ /= sumP_.last();
    sumP_ /= sumP_.last();
}


Foam::distributions::multiFixedValue::multiFixedValue
(
    const multiFixedValue& d,
    const label sampleQ
)
:
    FieldDistribution<distribution, multiFixedValue>(d, sampleQ),
    x_(d.x_),
    P_(d.P_),
    sumP_(d.sumP_)
{
    if (q() == d.q()) return;

    // Scale the probabilities
    P_ = integerPow(x_, q() - d.q())*d.P_;

    // Compute the cumulative sum of the probabilities
    sumP_.resize(x_.size() + 1);
    sumP_[0] = 0;
    forAll(x_, i)
    {
        sumP_[i + 1] = sumP_[i] + P_[i];
    }

    // Normalise
    P_ /= sumP_.last();
    sumP_ /= sumP_.last();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::distributions::multiFixedValue::~multiFixedValue()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::distributions::multiFixedValue::sample() const
{
    const scalar S = rndGen_.sample01<scalar>();

    label i = 0;
    for (; i < sumP_.size() - 1 && sumP_[i + 1] < S; ++ i);

    return x_[i];
}


Foam::scalar Foam::distributions::multiFixedValue::min() const
{
    return x_.first();
}


Foam::scalar Foam::distributions::multiFixedValue::max() const
{
    return x_.last();
}


Foam::scalar Foam::distributions::multiFixedValue::mean() const
{
    return sum(x_*P_);
}


Foam::tmp<Foam::scalarField>
Foam::distributions::multiFixedValue::CDF(const scalarField& x) const
{
    tmp<scalarField> tResult(new scalarField(x.size()));
    scalarField& result = tResult.ref();

    label i = 0;

    while (i < x.size() && x[i] < x_[0])
    {
        result[i] = 0;
        i ++;
    }

    for (label j = 0; j < x_.size() - 1; ++ j)
    {
        while (i < x.size() && x[i] < x_[j + 1])
        {
            result[i] = sumP_[j + 1];
            i ++;
        }
    }

    while (i < x.size())
    {
        result[i] = 1;
        i ++;
    }

    return tResult;
}


void Foam::distributions::multiFixedValue::write
(
    Ostream& os,
    const unitConversion& units
) const
{
    FieldDistribution<distribution, multiFixedValue>::write(os, units);

    // Recover the probabilities
    scalarField P(integerPow(x_, -q())*P_);

    // Normalise the probabilities
    P /= sum(P);

    // Construct and write the values
    List<Tuple2<scalar, scalar>> values(P_.size());
    forAll(values, i)
    {
        values[i].first() = units.toUser(x_[i]);
        values[i].second() = P[i];
    }
    writeEntry(os, "values", values);
}


Foam::tmp<Foam::scalarField>
Foam::distributions::multiFixedValue::plotX(const label) const
{
    tmp<scalarField> tResult(new scalarField(3*x_.size()));
    scalarField& result = tResult.ref();

    forAll(x_, i)
    {
        result[3*i] = result[3*i + 1] = result[3*i + 2] = x_[i];
    }

    return tResult;
}


Foam::tmp<Foam::scalarField>
Foam::distributions::multiFixedValue::plotPDF(const scalarField& x) const
{
    tmp<scalarField> tResult(new scalarField(3*x_.size()));
    scalarField& result = tResult.ref();

    forAll(x_, i)
    {
        result[3*i] = result[3*i + 2] = 0;
        result[3*i + 1] = P_[i]/rootVSmall;
    }

    return tResult;
}


// ************************************************************************* //
