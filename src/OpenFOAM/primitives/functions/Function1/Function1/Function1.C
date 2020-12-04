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

#include "Function1.H"


// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1<Type>::Function1(const word& name)
:
    name_(name)
{}


template<class Type>
Foam::Function1<Type>::Function1(const Function1<Type>& de)
:
    tmp<Function1<Type>>::refCount(),
    name_(de.name_)
{}


template<class Type, class Function1Type>
Foam::FieldFunction1<Type, Function1Type>::FieldFunction1
(
    const word& name
)
:
    Function1<Type>(name)
{}


template<class Type, class Function1Type>
Foam::tmp<Foam::Function1<Type>>
Foam::FieldFunction1<Type, Function1Type>::clone() const
{
    return tmp<Function1<Type>>
    (
        new Function1Type(refCast<const Function1Type>(*this))
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1<Type>::~Function1()
{}


template<class Type, class Function1Type>
Foam::FieldFunction1<Type, Function1Type>::~FieldFunction1()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Foam::word& Foam::Function1<Type>::name() const
{
    return name_;
}


template<class Type, class Function1Type>
Foam::tmp<Foam::Field<Type>> Foam::FieldFunction1<Type, Function1Type>::value
(
    const scalarField& x
) const
{
    tmp<Field<Type>> tfld(new Field<Type>(x.size()));
    Field<Type>& fld = tfld.ref();

    forAll(x, i)
    {
        fld[i] = refCast<const Function1Type>(*this).value(x[i]);
    }

    return tfld;
}


template<class Type, class Function1Type>
Foam::tmp<Foam::Field<Type>>
Foam::FieldFunction1<Type, Function1Type>::integral
(
    const scalarField& x1,
    const scalarField& x2
) const
{
    tmp<Field<Type>> tfld(new Field<Type>(x1.size()));
    Field<Type>& fld = tfld.ref();

    forAll(x1, i)
    {
        fld[i] = refCast<const Function1Type>(*this).integral(x1[i], x2[i]);
    }

    return tfld;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::Function1<Type>::operator=(const Function1<Type>& f)
{
    if (this == &f)
    {
        FatalErrorInFunction
            << "attempted assignment to self"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * IOstream Functions  * * * * * * * * * * * * //

template<class Type>
void  Foam::writeEntry(Ostream& os, const Function1<Type>& f1)
{
    writeKeyword(os, f1.name())
        << nl << indent << token::BEGIN_BLOCK << nl << incrIndent;

    writeEntry(os, "type", f1.type());

    f1.write(os);

    os  << decrIndent << indent << token::END_BLOCK << endl;
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

    f1.write(os);

    return os;
}


// ************************************************************************* //
