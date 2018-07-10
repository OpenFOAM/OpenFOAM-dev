/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "TableBase.H"
#include "Time.H"
#include "interpolationWeights.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
const Foam::interpolationWeights&
Foam::Function1Types::TableBase<Type>::interpolator() const
{
    if (interpolatorPtr_.empty())
    {
        // Re-work table into linear list
        tableSamplesPtr_.reset(new scalarField(table_.size()));
        scalarField& tableSamples = tableSamplesPtr_();
        forAll(table_, i)
        {
            tableSamples[i] = table_[i].first();
        }
        interpolatorPtr_ = interpolationWeights::New
        (
            interpolationScheme_,
            tableSamples
        );
    }

    return interpolatorPtr_();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::TableBase<Type>::TableBase
(
    const word& name,
    const dictionary& dict
)
:
    Function1<Type>(name),
    name_(name),
    boundsHandling_
    (
        wordToBoundsHandling
        (
            dict.lookupOrDefault<word>("outOfBounds", "clamp")
        )
    ),
    interpolationScheme_
    (
        dict.lookupOrDefault<word>("interpolationScheme", "linear")
    ),
    table_()
{}


template<class Type>
Foam::Function1Types::TableBase<Type>::TableBase(const TableBase<Type>& tbl)
:
    Function1<Type>(tbl),
    name_(tbl.name_),
    boundsHandling_(tbl.boundsHandling_),
    interpolationScheme_(tbl.interpolationScheme_),
    table_(tbl.table_),
    tableSamplesPtr_(tbl.tableSamplesPtr_),
    interpolatorPtr_(tbl.interpolatorPtr_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::TableBase<Type>::~TableBase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::word Foam::Function1Types::TableBase<Type>::boundsHandlingToWord
(
     const boundsHandling& bound
) const
{
    word enumName("warn");

    switch (bound)
    {
        case ERROR:
        {
            enumName = "error";
            break;
        }
        case WARN:
        {
            enumName = "warn";
            break;
        }
        case CLAMP:
        {
            enumName = "clamp";
            break;
        }
        case REPEAT:
        {
            enumName = "repeat";
            break;
        }
    }

    return enumName;
}


template<class Type>
typename Foam::Function1Types::TableBase<Type>::boundsHandling
Foam::Function1Types::TableBase<Type>::wordToBoundsHandling
(
    const word& bound
) const
{
    if (bound == "error")
    {
        return ERROR;
    }
    else if (bound == "warn")
    {
        return WARN;
    }
    else if (bound == "clamp")
    {
        return CLAMP;
    }
    else if (bound == "repeat")
    {
        return REPEAT;
    }
    else
    {
        WarningInFunction
            << "bad outOfBounds specifier " << bound << " using 'warn'"
            << endl;

        return WARN;
    }
}


template<class Type>
typename Foam::Function1Types::TableBase<Type>::boundsHandling
Foam::Function1Types::TableBase<Type>::outOfBounds
(
    const boundsHandling& bound
)
{
    boundsHandling prev = boundsHandling_;
    boundsHandling_ = bound;

    return prev;
}


template<class Type>
void Foam::Function1Types::TableBase<Type>::check() const
{
    if (!table_.size())
    {
        FatalErrorInFunction
            << "Table for entry " << this->name_ << " is invalid (empty)"
            << nl << exit(FatalError);
    }

    label n = table_.size();
    scalar prevValue = table_[0].first();

    for (label i = 1; i < n; ++i)
    {
        const scalar currValue = table_[i].first();

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
bool Foam::Function1Types::TableBase<Type>::checkMinBounds
(
    const scalar x,
    scalar& xDash
) const
{
    if (x < table_[0].first())
    {
        switch (boundsHandling_)
        {
            case ERROR:
            {
                FatalErrorInFunction
                    << "value (" << x << ") underflow"
                    << exit(FatalError);
                break;
            }
            case WARN:
            {
                WarningInFunction
                    << "value (" << x << ") underflow" << nl
                    << endl;

                // fall-through to 'CLAMP'
                [[fallthrough]];
            }
            case CLAMP:
            {
                xDash = table_[0].first();
                return true;
                break;
            }
            case REPEAT:
            {
                // adjust x to >= minX
                scalar span = table_.last().first() - table_[0].first();
                xDash = fmod(x - table_[0].first(), span) + table_[0].first();
                break;
            }
        }
    }
    else
    {
        xDash = x;
    }

    return false;
}


template<class Type>
bool Foam::Function1Types::TableBase<Type>::checkMaxBounds
(
    const scalar x,
    scalar& xDash
) const
{
    if (x > table_.last().first())
    {
        switch (boundsHandling_)
        {
            case ERROR:
            {
                FatalErrorInFunction
                    << "value (" << x << ") overflow"
                    << exit(FatalError);
                break;
            }
            case WARN:
            {
                WarningInFunction
                    << "value (" << x << ") overflow" << nl
                    << endl;

                // fall-through to 'CLAMP'
                [[fallthrough]];
            }
            case CLAMP:
            {
                xDash = table_.last().first();
                return true;
                break;
            }
            case REPEAT:
            {
                // adjust x to >= minX
                scalar span = table_.last().first() - table_[0].first();
                xDash = fmod(x - table_[0].first(), span) + table_[0].first();
                break;
            }
        }
    }
    else
    {
        xDash = x;
    }

    return false;
}


template<class Type>
void Foam::Function1Types::TableBase<Type>::convertTimeBase(const Time& t)
{
    forAll(table_, i)
    {
        scalar value = table_[i].first();
        table_[i].first() = t.userTimeToTime(value);
    }

    tableSamplesPtr_.clear();
    interpolatorPtr_.clear();
}


template<class Type>
Type Foam::Function1Types::TableBase<Type>::value(const scalar x) const
{
    scalar xDash = x;

    if (checkMinBounds(x, xDash))
    {
        return table_[0].second();
    }

    if (checkMaxBounds(xDash, xDash))
    {
        return table_.last().second();
    }

    // Use interpolator
    interpolator().valueWeights(xDash, currentIndices_, currentWeights_);

    Type t = currentWeights_[0]*table_[currentIndices_[0]].second();
    for (label i = 1; i < currentIndices_.size(); i++)
    {
        t += currentWeights_[i]*table_[currentIndices_[i]].second();
    }

    return t;
}


template<class Type>
Type Foam::Function1Types::TableBase<Type>::integrate
(
    const scalar x1,
    const scalar x2
) const
{
    // Use interpolator
    interpolator().integrationWeights(x1, x2, currentIndices_, currentWeights_);

    Type sum = currentWeights_[0]*table_[currentIndices_[0]].second();
    for (label i = 1; i < currentIndices_.size(); i++)
    {
       sum += currentWeights_[i]*table_[currentIndices_[i]].second();
    }

    return sum;
}


template<class Type>
Foam::tmp<Foam::scalarField> Foam::Function1Types::TableBase<Type>::x() const
{
    tmp<scalarField> tfld(new scalarField(table_.size(), 0.0));
    scalarField& fld = tfld.ref();

    forAll(table_, i)
    {
        fld[i] = table_[i].first();
    }

    return tfld;
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::Function1Types::TableBase<Type>::y() const
{
    tmp<Field<Type>> tfld(new Field<Type>(table_.size(), Zero));
    Field<Type>& fld = tfld.ref();

    forAll(table_, i)
    {
        fld[i] = table_[i].second();
    }

    return tfld;
}


template<class Type>
void Foam::Function1Types::TableBase<Type>::writeEntries(Ostream& os) const
{
    if (boundsHandling_ != CLAMP)
    {
        os.writeKeyword("outOfBounds") << boundsHandlingToWord(boundsHandling_)
            << token::END_STATEMENT << nl;
    }
    if (interpolationScheme_ != "linear")
    {
        os.writeKeyword("interpolationScheme") << interpolationScheme_
            << token::END_STATEMENT << nl;
    }
}


template<class Type>
void Foam::Function1Types::TableBase<Type>::writeData(Ostream& os) const
{
    Function1<Type>::writeData(os);
    os  << nl << indent << table_ << token::END_STATEMENT << nl;
    writeEntries(os);
}


// ************************************************************************* //
