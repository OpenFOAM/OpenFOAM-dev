/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
Foam::scalar Foam::Function1Types::TableBase<Type>::bound(const scalar x) const
{
    const bool under = x < table_.first().first();
    const bool over = x > table_.last().first();

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
                const scalar t0 = table_.first().first();
                const scalar t1 = table_.last().first();
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
Foam::Function1Types::TableBase<Type>::TableBase
(
    const word& name,
    const dictionary& dict
)
:
    tableBase(),
    Function1<Type>(name),
    name_(name),
    boundsHandling_
    (
        dict.found("outOfBounds")
      ? tableBase::boundsHandlingNames_.read(dict.lookup("outOfBounds"))
      : tableBase::boundsHandling::warn
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
    tableBase(),
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
Type Foam::Function1Types::TableBase<Type>::value(const scalar x) const
{
    const scalar bx = bound(x);

    Type y = Zero;

    interpolator().valueWeights(bx, indices_, weights_);
    forAll(indices_, i)
    {
        y += weights_[i]*table_[indices_[i]].second();
    }

    return y;
}


template<class Type>
Type Foam::Function1Types::TableBase<Type>::integrate
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
       sumY += weights_[i]*table_[indices_[i]].second();
    }

    if (boundsHandling_ == tableBase::boundsHandling::repeat)
    {
        const scalar t0 = table_.first().first();
        const scalar t1 = table_.last().first();
        const scalar dt = t1 - t0;
        const label n = floor((x2 - t0)/dt) - floor((x1 - t0)/dt);

        if (n != 0)
        {
            Type sumY01 = Zero;

            interpolator().integrationWeights(t0, t1, indices_, weights_);

            forAll(indices_, i)
            {
                sumY01 += weights_[i]*table_[indices_[i]].second();
            }
            sumY += n*sumY01;
        }
    }

    return sumY;
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
    if (boundsHandling_ != tableBase::boundsHandling::clamp)
    {
        writeEntry
        (
            os,
            "outOfBounds",
            tableBase::boundsHandlingNames_[boundsHandling_]
        );
    }
    if (interpolationScheme_ != "linear")
    {
        writeEntry(os, "interpolationScheme", interpolationScheme_);
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
