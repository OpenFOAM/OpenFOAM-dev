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


template<class Function1Type>
Foam::tmp<Foam::Field<typename Function1Type::returnType>>
Foam::FieldFunction1<Function1Type>::value
(
    const scalarField& x
) const
{
    tmp<Field<Type>> tfld(new Field<Type>(x.size()));
    Field<Type>& fld = tfld.ref();

    forAll(x, i)
    {
        fld[i] = Function1Type::value(x[i]);
    }
    return tfld;
}


template<class Function1Type>
Foam::FieldFunction1<Function1Type>::FieldFunction1
(
    const word& entryName,
    const dictionary& dict
)
:
    Function1Type(entryName, dict)
{}


template<class Function1Type>
Foam::tmp<Foam::Function1<typename Function1Type::returnType>>
Foam::FieldFunction1<Function1Type>::clone() const
{
    return tmp<Function1<Type>>
    (
        new FieldFunction1<Function1Type>(*this)
    );
}


template<class Function1Type>
Foam::tmp<Foam::Field<typename Function1Type::returnType>>
Foam::FieldFunction1<Function1Type>::integrate
(
    const scalarField& x1,
    const scalarField& x2
) const
{
    tmp<Field<Type>> tfld(new Field<Type>(x1.size()));
    Field<Type>& fld = tfld.ref();

    forAll(x1, i)
    {
        fld[i] = Function1Type::integrate(x1[i], x2[i]);
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
