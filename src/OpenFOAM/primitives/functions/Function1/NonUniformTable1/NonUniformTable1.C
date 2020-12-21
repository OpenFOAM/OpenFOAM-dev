/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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

#include "NonUniformTable1.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1s::NonUniformTable<Type>::NonUniformTable
(
    const word& name,
    const dictionary& dict
)
:
    FieldFunction1<Type, NonUniformTable<Type>>(name),
    low_(great),
    high_(-great),
    values_(),
    delta_(great),
    reader_(TableReader<Type>::New(name, dict, this->values_))
{
    if (values_.size() < 2)
    {
        FatalIOErrorInFunction(dict)
            << "Table " << nl
            << "    " << name << nl
            << "    has less than 2 entries."
            << exit(FatalIOError);
    }
    else
    {
        low_ = values_.first().first();
        high_ = values_.last().first();

        for(label i = 1; i<values_.size(); i++)
        {
            delta_ = min(delta_, values_[i].first() - values_[i - 1].first());
        }

        delta_ *= 0.9;

        jumpTable_.setSize((high_ - low_)/delta_ + 1);

        label i = 0;
        forAll(jumpTable_, j)
        {
            const scalar x = low_ + j*delta_;

            if (x > values_[i + 1].first())
            {
                i++;
            }

            jumpTable_[j] = i;
        }
    }
}


template<class Type>
Foam::Function1s::NonUniformTable<Type>::NonUniformTable
(
    const NonUniformTable<Type>& nut
)
:
    FieldFunction1<Type, NonUniformTable<Type>>(nut),
    low_(nut.low_),
    high_(nut.high_),
    values_(nut.values_),
    delta_(nut.delta_),
    jumpTable_(nut.jumpTable_),
    reader_(nut.reader_, false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::Function1s::NonUniformTable<Type>::value
(
    scalar x
) const
{
    const label i = index(x);
    const scalar xi = values_[i].first();
    const scalar lambda = (x - xi)/(values_[i + 1].first() - xi);

    return
        values_[i].second()
      + lambda*(values_[i + 1].second() - values_[i].second());
}


template<class Type>
Type Foam::Function1s::NonUniformTable<Type>::integral
(
    const scalar x1,
    const scalar x2
) const
{
    NotImplemented;
    return Zero;
}


template<class Type>
Type Foam::Function1s::NonUniformTable<Type>::dfdT
(
    scalar T
) const
{
    const label i = index(T);

    return
        (values_[i + 1].second() - values_[i].second())
       /(values_[i + 1].first() - values_[i].first());
}


template<class Type>
void Foam::Function1s::NonUniformTable<Type>::write(Ostream& os) const
{
    reader_->write(os, values_);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::Function1s::NonUniformTable<Type>::operator=
(
    const NonUniformTable<Type>& nut
)
{
    low_ = nut.low_;
    high_ = nut.high_;
    values_ = nut.values_;
    delta_ = nut.delta_;
    jumpTable_ = nut.jumpTable_;
    reader_ = nut.reader_->clone();
}


// ************************************************************************* //
