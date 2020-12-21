/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
Foam::scalar Foam::Function1s::Table<Type>::bound
(
    const scalar x
) const
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
            case tableBase::boundsHandling::clamp:
            {
                break;
            }
            case tableBase::boundsHandling::repeat:
            {
                const scalar t0 = values_.first().first();
                const scalar t1 = values_.last().first();
                const scalar dt = t1 - t0;
                const label n = floor((x - t0)/dt);
                return x - n*dt;
            }
        }
    }

    return x;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1s::Table<Type>::Table
(
    const word& name,
    const tableBase::boundsHandling boundsHandling,
    const word& interpolationScheme,
    const List<Tuple2<scalar, Type>>& table
)
:
    FieldFunction1<Type, Table<Type>>(name),
    boundsHandling_(boundsHandling),
    interpolationScheme_(interpolationScheme),
    values_(table)
{}


template<class Type>
Foam::Function1s::Table<Type>::Table
(
    const word& name,
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
    values_(),
    reader_(TableReader<Type>::New(name, dict, this->values_))
{
    check();
}


template<class Type>
Foam::Function1s::Table<Type>::Table(const Table<Type>& tbl)
:
    FieldFunction1<Type, Table<Type>>(tbl),
    boundsHandling_(tbl.boundsHandling_),
    interpolationScheme_(tbl.interpolationScheme_),
    values_(tbl.values_),
    tableSamplesPtr_(tbl.tableSamplesPtr_),
    interpolatorPtr_(tbl.interpolatorPtr_),
    reader_(tbl.reader_, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1s::Table<Type>::~Table()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::Function1s::Table<Type>::value
(
    const scalar x
) const
{
    const scalar bx = bound(x);

    Type y = Zero;

    interpolator().valueWeights(bx, indices_, weights_);
    forAll(indices_, i)
    {
        y += weights_[i]*values_[indices_[i]].second();
    }

    return y;
}


template<class Type>
Type Foam::Function1s::Table<Type>::integral
(
    const scalar x1,
    const scalar x2
) const
{
    const scalar bx1 = bound(x1), bx2 = bound(x2);

    Type sumY = Zero;

    interpolator().integrationWeights(bx1, bx2, indices_, weights_);
    forAll(indices_, i)
    {
       sumY += weights_[i]*values_[indices_[i]].second();
    }

    if (boundsHandling_ == tableBase::boundsHandling::repeat)
    {
        const scalar t0 = values_.first().first();
        const scalar t1 = values_.last().first();
        const scalar dt = t1 - t0;
        const label n = floor((x2 - t0)/dt) - floor((x1 - t0)/dt);

        if (n != 0)
        {
            Type sumY01 = Zero;

            interpolator().integrationWeights(t0, t1, indices_, weights_);

            forAll(indices_, i)
            {
                sumY01 += weights_[i]*values_[indices_[i]].second();
            }
            sumY += n*sumY01;
        }
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
void Foam::Function1s::Table<Type>::write(Ostream& os) const
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

    reader_->write(os, values_);
}


// ************************************************************************* //
