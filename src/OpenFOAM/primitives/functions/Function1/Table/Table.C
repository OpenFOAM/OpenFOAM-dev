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

#include "Table.H"
#include "linearInterpolationWeights.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
const Foam::interpolationWeights&
Foam::Function1s::Table<Type>::interpolator() const
{
    if (interpolatorPtr_.empty())
    {
        // Re-work table into linear list
        tableSamplesPtr_.reset(new scalarField(values_.size()));
        scalarField& tableSamples = tableSamplesPtr_();
        forAll(values_, i)
        {
            tableSamples[i] = values_[i].first();
        }
        interpolatorPtr_ = interpolationWeights::New
        (
            interpolationScheme_,
            tableSamples
        );
    }

    return interpolatorPtr_();
}


template<class Type>
void Foam::Function1s::Table<Type>::check() const
{
    if (!values_.size())
    {
        FatalErrorInFunction
            << "Table for entry " << this->name() << " is invalid (empty)"
            << nl << exit(FatalError);
    }

    label n = values_.size();
    scalar prevValue = values_[0].first();

    for (label i = 1; i < n; ++i)
    {
        const scalar currValue = values_[i].first();

        // avoid duplicate values (divide-by-zero error)
        if (currValue <= prevValue)
        {
            FatalErrorInFunction
                << "out-of-order value: " << currValue << " at index " << i
                << exit(FatalError);
        }
        prevValue = currValue;
    }
}


template<class Type>
void Foam::Function1s::Table<Type>::checkX(const scalar x) const
{
    const bool under = x < values_.first().first();
    const bool over = x > values_.last().first();

    auto errorMessage = [&]()
    {
        return "value (" + name(x) + ") " + (under ? "under" : "over") + "flow";
    };

    if (under || over)
    {
        switch (boundsHandling_)
        {
            case tableBase::boundsHandling::error:
            {
                FatalErrorInFunction
                    << errorMessage() << nl << exit(FatalError);
                break;
            }
            case tableBase::boundsHandling::warn:
            {
                WarningInFunction
                    << errorMessage() << nl << endl;
                break;
            }
            default:
            {
                break;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1s::Table<Type>::Table
(
    const word& name,
    const tableBase::boundsHandling boundsHandling,
    const word& interpolationScheme,
    const autoPtr<TableReader<Type>>& reader,
    const List<Tuple2<scalar, Type>>& table
)
:
    FieldFunction1<Type, Table<Type>>(name),
    boundsHandling_(boundsHandling),
    interpolationScheme_(interpolationScheme),
    reader_(reader, false),
    values_(table)
{}


template<class Type>
Foam::Function1s::Table<Type>::Table
(
    const word& name,
    const unitConversions& units,
    const dictionary& dict
)
:
    FieldFunction1<Type, Table<Type>>(name),
    boundsHandling_
    (
        dict.found("outOfBounds")
      ? tableBase::boundsHandlingNames_.read(dict.lookup("outOfBounds"))
      : tableBase::boundsHandling::clamp
    ),
    interpolationScheme_
    (
        dict.lookupOrDefault<word>
        (
            "interpolationScheme",
            linearInterpolationWeights::typeName
        )
    ),
    reader_(TableReader<Type>::New(name, units, dict)),
    values_(reader_->read(units, dict))
{
    check();
}


template<class Type>
Foam::Function1s::Table<Type>::Table
(
    const word& name,
    const unitConversion& xUnits,
    const unitConversion& valueUnits,
    const dictionary& dict
)
:
    Table(name, {xUnits, valueUnits}, dict)
{}


template<class Type>
Foam::Function1s::Table<Type>::Table
(
    const word& name,
    const unitConversions& units,
    Istream& is
)
:
    FieldFunction1<Type, Table<Type>>(name),
    boundsHandling_(tableBase::boundsHandling::clamp),
    interpolationScheme_(linearInterpolationWeights::typeName),
    reader_(new TableReaders::Embedded<Type>()),
    values_(TableReaders::Embedded<Type>().read(units, is))
{
    check();
}


template<class Type>
Foam::Function1s::Table<Type>::Table(const Table<Type>& tbl)
:
    FieldFunction1<Type, Table<Type>>(tbl),
    boundsHandling_(tbl.boundsHandling_),
    interpolationScheme_(tbl.interpolationScheme_),
    reader_(tbl.reader_, false),
    values_(tbl.values_),
    tableSamplesPtr_(tbl.tableSamplesPtr_),
    interpolatorPtr_(tbl.interpolatorPtr_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1s::Table<Type>::~Table()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::Function1s::Table<Type>::value
(
    const scalar xArg
) const
{
    scalar x(xArg);

    checkX(x);

    const scalar x0 = values_.first().first();
    const scalar x1 = values_.last().first();

    Type y = Zero;

    switch (boundsHandling_)
    {
        case tableBase::boundsHandling::clamp:
        {
            // Integration weights clamp by default, so do nothing
            break;
        }
        case tableBase::boundsHandling::zero:
        {
            // Return zero if outside the range
            if (x <= x0 || x >= x1)
            {
                return Zero;
            }
            break;
        }
        case tableBase::boundsHandling::repeat:
        {
            // Shift the value into the range of the table
            const scalar n = floor((x - x0)/(x1 - x0));
            x -= n*(x1 - x0);
            break;
        }
        default:
        {
            break;
        }
    }

    // Evaluate
    interpolator().valueWeights(x, indices_, weights_);
    forAll(indices_, i)
    {
        y += weights_[i]*values_[indices_[i]].second();
    }

    return y;
}


template<class Type>
Type Foam::Function1s::Table<Type>::integral
(
    const scalar xAarg,
    const scalar xBarg
) const
{
    scalar xA(xAarg), xB(xBarg);

    checkX(xA);
    checkX(xB);

    const scalar x0 = values_.first().first();
    const scalar x1 = values_.last().first();

    Type sumY = Zero;

    switch (boundsHandling_)
    {
        case tableBase::boundsHandling::clamp:
        {
            // Integration weights clamp by default, so do nothing
            break;
        }
        case tableBase::boundsHandling::zero:
        {
            // Actually clamp to remove any integral components from the ends
            xA = min(max(xA, x0), x1);
            xB = min(max(xB, x0), x1);
            break;
        }
        case tableBase::boundsHandling::repeat:
        {
            // Shift the values into the range of the table and add any
            // repetitions necessary to the integral
            const scalar nA = floor((xA - x0)/(x1 - x0));
            const scalar nB = floor((xB - x0)/(x1 - x0));
            if (nA != nB)
            {
                interpolator().integrationWeights(x0, x1, indices_, weights_);
                forAll(indices_, i)
                {
                    sumY += (nB - nA)*weights_[i]*values_[indices_[i]].second();
                }
            }
            xA -= nA*(x1 - x0);
            xB -= nB*(x1 - x0);
            break;
        }
        default:
        {
            break;
        }
    }

    // Integrate between the bounds
    interpolator().integrationWeights(xA, xB, indices_, weights_);
    forAll(indices_, i)
    {
       sumY += weights_[i]*values_[indices_[i]].second();
    }

    return sumY;
}


template<class Type>
Foam::tmp<Foam::scalarField>
Foam::Function1s::Table<Type>::x() const
{
    tmp<scalarField> tfld(new scalarField(values_.size(), 0.0));
    scalarField& fld = tfld.ref();

    forAll(values_, i)
    {
        fld[i] = values_[i].first();
    }

    return tfld;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::Function1s::Table<Type>::y() const
{
    tmp<Field<Type>> tfld(new Field<Type>(values_.size(), Zero));
    Field<Type>& fld = tfld.ref();

    forAll(values_, i)
    {
        fld[i] = values_[i].second();
    }

    return tfld;
}


template<class Type>
void Foam::Function1s::Table<Type>::write
(
    Ostream& os,
    const unitConversions& units
) const
{
    writeEntryIfDifferent
    (
        os,
        "outOfBounds",
        tableBase::boundsHandlingNames_[tableBase::boundsHandling::clamp],
        tableBase::boundsHandlingNames_[boundsHandling_]
    );

    writeEntryIfDifferent
    (
        os,
        "interpolationScheme",
        linearInterpolationWeights::typeName,
        interpolationScheme_
    );

    reader_->write(os, units, values_);
}


// ************************************************************************* //
