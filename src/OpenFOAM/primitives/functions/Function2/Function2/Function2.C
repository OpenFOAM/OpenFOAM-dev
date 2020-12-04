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

#include "Function2.H"

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function2<Type>::Function2(const word& name)
:
    name_(name)
{}


template<class Type>
Foam::Function2<Type>::Function2(const Function2<Type>& de)
:
    tmp<Function2<Type>>::refCount(),
    name_(de.name_)
{}


template<class Type, class Function2Type>
Foam::FieldFunction2<Type, Function2Type>::FieldFunction2
(
    const word& name
)
:
    Function2<Type>(name)
{}


template<class Type, class Function2Type>
Foam::tmp<Foam::Function2<Type>>
Foam::FieldFunction2<Type, Function2Type>::clone() const
{
    return tmp<Function2<Type>>
    (
        new Function2Type(refCast<const Function2Type>(*this))
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function2<Type>::~Function2()
{}


template<class Type, class Function2Type>
Foam::FieldFunction2<Type, Function2Type>::~FieldFunction2()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Foam::word& Foam::Function2<Type>::name() const
{
    return name_;
}


template<class Type, class Function2Type>
Foam::tmp<Foam::Field<Type>> Foam::FieldFunction2<Type, Function2Type>::value
(
    const scalarField& x,
    const scalarField& y
) const
{
    tmp<Field<Type>> tfld(new Field<Type>(x.size()));
    Field<Type>& fld = tfld.ref();

    forAll(x, i)
    {
        fld[i] = refCast<const Function2Type>(*this).value(x[i], y[i]);
    }

    return tfld;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::Function2<Type>::operator=(const Function2<Type>& f)
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
void  Foam::writeEntry(Ostream& os, const Function2<Type>& f2)
{
    writeKeyword(os, f2.name())
        << nl << indent << token::BEGIN_BLOCK << nl << incrIndent;

    writeEntry(os, "type", f2.type());

    f2.write(os);

    os  << decrIndent << indent << token::END_BLOCK << endl;
}


// * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const Function2<Type>& f2
)
{
    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const Function2<Type>&)"
    );

    f2.write(os);

    return os;
}


// ************************************************************************* //
