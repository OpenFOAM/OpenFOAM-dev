/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "DynamicField.H"
#include "expressionAssert.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
Foam::DynamicField<T, SizeInc, SizeMult, SizeDiv>::DynamicField(Istream& is)
:
    Field<T>(is),
    capacity_(Field<T>::size())
{}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
Foam::DynamicField<T, SizeInc, SizeMult, SizeDiv>::DynamicField
(
    const word& keyword,
    const dictionary& dict,
    const label size
)
:
    Field<T>(keyword, dict, size),
    capacity_(Field<T>::size())
{}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
Foam::DynamicField<T, SizeInc, SizeMult, SizeDiv>::DynamicField
(
    const word& keyword,
    const unitSet& defaultUnits,
    const dictionary& dict,
    const label size
)
:
    Field<T>(keyword, defaultUnits, dict, size),
    capacity_(Field<T>::size())
{}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
Foam::tmp<Foam::DynamicField<T, SizeInc, SizeMult, SizeDiv>>
Foam::DynamicField<T, SizeInc, SizeMult, SizeDiv>::clone() const
{
    return tmp<DynamicField<T, SizeInc, SizeMult, SizeDiv>>
    (
        new DynamicField<T, SizeInc, SizeMult, SizeDiv>(*this)
    );
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
template<class Expression, class>
void Foam::DynamicField<T, SizeInc, SizeMult, SizeDiv>::operator=
(
    const Expression& e
)
{
    #ifdef FULLDEBUG
    expression::assertSameAllContainerProperty<expression::Size>(e);
    #endif

    // Get size of the first field in the expression
    const label size = expression::getAll<expression::Size>(e);

    if (capacity_ >= size)
    {
        // Can copy w/o reallocating, match initial size to avoid reallocation
        Field<T>::size(size);
        Field<T>::operator=(e);
    }
    else
    {
        // Make everything available for the copy operation
        Field<T>::size(capacity_);

        Field<T>::operator=(e);
        capacity_ = Field<T>::size();
    }
}


// * * * * * * * * * * * * * * * IOstream Operator * * * * * * * * * * * * * //

template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const DynamicField<T, SizeInc, SizeMult, SizeDiv>& lst
)
{
    os << static_cast<const Field<T>&>(lst);
    return os;
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
Foam::Istream& Foam::operator>>
(
    Istream& is,
    DynamicField<T, SizeInc, SizeMult, SizeDiv>& lst
)
{
    is >> static_cast<Field<T>&>(lst);
    lst.capacity_ = lst.Field<T>::size();

    return is;
}


// ************************************************************************* //
