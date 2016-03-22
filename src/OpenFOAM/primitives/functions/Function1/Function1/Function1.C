/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "Function1.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1<Type>::Function1(const word& entryName)
:
    name_(entryName)
{}


template<class Type>
Foam::Function1<Type>::Function1(const Function1<Type>& de)
:
    tmp<Function1<Type>>::refCount(),
    name_(de.name_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1<Type>::~Function1()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Foam::word& Foam::Function1<Type>::name() const
{
    return name_;
}


template<class Type>
void Foam::Function1<Type>::convertTimeBase(const Time&)
{}


template<class Type>
Type Foam::Function1<Type>::value(const scalar x) const
{
    NotImplemented;

    return Zero;
}


template<class Type>
Type Foam::Function1<Type>::integrate(const scalar x1, const scalar x2) const
{
    NotImplemented;

    return Zero;
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::Function1<Type>::value
(
    const scalarField& x
) const
{
    tmp<Field<Type>> tfld(new Field<Type>(x.size()));
    Field<Type>& fld = tfld.ref();

    forAll(x, i)
    {
        fld[i] = this->value(x[i]);
    }
    return tfld;
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::Function1<Type>::integrate
(
    const scalarField& x1,
    const scalarField& x2
) const
{
    tmp<Field<Type>> tfld(new Field<Type>(x1.size()));
    Field<Type>& fld = tfld.ref();

    forAll(x1, i)
    {
        fld[i] = this->integrate(x1[i], x2[i]);
    }
    return tfld;
}


template<class Type>
void Foam::Function1<Type>::writeData(Ostream& os) const
{
    os.writeKeyword(name_) << type();
}


// * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const Function1<Type>& f1
)
{
    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const Function1<Type>&)"
    );

    os  << f1.name_;
    f1.writeData(os);

    return os;
}


// ************************************************************************* //
