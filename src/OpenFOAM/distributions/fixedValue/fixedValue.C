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

#include "fixedValue.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace distributions
{
    defineTypeNameAndDebug(fixedValue, 0);
    addToRunTimeSelectionTable(distribution, fixedValue, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distributions::fixedValue::fixedValue
(
    const unitConversion& units,
    const dictionary& dict,
    const label,
    randomGenerator&& rndGen
)
:
    FieldDistribution<distribution, fixedValue>
    (
        -labelMax,
        -labelMax,
        std::move(rndGen)
    ),
    value_(dict.lookup<scalar>("value", units))
{}


Foam::distributions::fixedValue::fixedValue
(
    const fixedValue& d,
    const label sampleQ
)
:
    FieldDistribution<distribution, fixedValue>(d, sampleQ),
    value_(d.value_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::distributions::fixedValue::~fixedValue()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::distributions::fixedValue::sample() const
{
    return value_;
}


Foam::scalar Foam::distributions::fixedValue::min() const
{
    return value_;
}


Foam::scalar Foam::distributions::fixedValue::max() const
{
    return value_;
}


Foam::scalar Foam::distributions::fixedValue::mean() const
{
    return value_;
}


Foam::tmp<Foam::scalarField>
Foam::distributions::fixedValue::CDF(const scalarField& x) const
{
    return pos(x - value_);
}


void Foam::distributions::fixedValue::write
(
    Ostream& os,
    const unitConversion& units
) const
{
    writeEntry(os, "type", type());
    writeEntry(os, "value", units, value_);
}


Foam::tmp<Foam::scalarField>
Foam::distributions::fixedValue::plotX(const label n) const
{
    const scalar d = 0.1*mag(value_);

    tmp<scalarField> tResult(new scalarField(5));
    scalarField& result = tResult.ref();

    result[0] = value_ - d;
    result[1] = value_*(1 - sign(value_)*rootVSmall);
    result[2] = value_;
    result[3] = value_*(1 + sign(value_)*rootVSmall);
    result[4] = value_ + d;

    return tResult;
}


Foam::tmp<Foam::scalarField>
Foam::distributions::fixedValue::plotPDF(const scalarField& x) const
{
    tmp<scalarField> tResult(new scalarField(5, 0));
    scalarField& result = tResult.ref();

    result[2] = 1/rootVSmall;

    return tResult;
}


// ************************************************************************* //
