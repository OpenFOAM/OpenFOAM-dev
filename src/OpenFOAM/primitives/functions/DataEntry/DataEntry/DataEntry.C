/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "DataEntry.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

template<class Type>
Foam::DataEntry<Type>::DataEntry(const word& entryName)
:
    refCount(),
    name_(entryName)
{}


template<class Type>
Foam::DataEntry<Type>::DataEntry(const DataEntry<Type>& de)
:
    refCount(),
    name_(de.name_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::DataEntry<Type>::~DataEntry()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Foam::word& Foam::DataEntry<Type>::name() const
{
    return name_;
}


template<class Type>
void Foam::DataEntry<Type>::convertTimeBase(const Time&)
{
    // do nothing
}


template<class Type>
Type Foam::DataEntry<Type>::value(const scalar x) const
{
    notImplemented("Type Foam::DataEntry<Type>::value(const scalar) const");

    return pTraits<Type>::zero;
}


template<class Type>
Type Foam::DataEntry<Type>::integrate(const scalar x1, const scalar x2) const
{
    notImplemented
    (
        "Type Foam::DataEntry<Type>::integrate"
        "("
            "const scalar, "
            "const scalar"
        ") const"
    );

    return pTraits<Type>::zero;
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::DataEntry<Type>::value
(
    const scalarField& x
) const
{
    tmp<Field<Type> > tfld(new Field<Type>(x.size()));
    Field<Type>& fld = tfld();

    forAll(x, i)
    {
        fld[i] = this->value(x[i]);
    }
    return tfld;
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::DataEntry<Type>::integrate
(
    const scalarField& x1,
    const scalarField& x2
) const
{
    tmp<Field<Type> > tfld(new Field<Type>(x1.size()));
    Field<Type>& fld = tfld();

    forAll(x1, i)
    {
        fld[i] = this->integrate(x1[i], x2[i]);
    }
    return tfld;
}



template<class Type>
Foam::dimensioned<Type> Foam::DataEntry<Type>::dimValue(const scalar x) const
{
    notImplemented
    (
        "dimensioned<Type> Foam::DataEntry<dimensioned<Type> >::dimValue"
        "(const scalar) const"
    );

    return dimensioned<Type>("zero", dimless, pTraits<Type>::zero);
}


template<class Type>
Foam::dimensioned<Type> Foam::DataEntry<Type>::dimIntegrate
(
    const scalar x1,
    const scalar x2
) const
{
    notImplemented
    (
        "dimensioned<Type> Foam::DataEntry<Type>::dimIntegrate"
        "("
            "const scalar, "
            "const scalar"
        ") const"
    );

    return dimensioned<Type>("zero", dimless, pTraits<Type>::zero);
}


template<class Type>
Foam::tmp<Foam::Field<Foam::dimensioned<Type> > >
Foam::DataEntry<Type>::dimValue
(
    const scalarField& x
) const
{

    tmp<Field<dimensioned<Type> > > tfld
    (
        new Field<dimensioned<Type> >
        (
            x.size(),
            dimensioned<Type>("zero", dimless, pTraits<Type>::zero)
        )
    );

    Field<dimensioned<Type> >& fld = tfld();

    forAll(x, i)
    {
        fld[i] = this->dimValue(x[i]);
    }
    return tfld;
}


template<class Type>
Foam::tmp<Foam::Field<Foam::dimensioned<Type> > >
Foam::DataEntry<Type>::dimIntegrate
(
    const scalarField& x1,
    const scalarField& x2
) const
{
    tmp<Field<dimensioned<Type> > > tfld
    (
        new Field<dimensioned<Type> >(x1.size())
    );

    Field<dimensioned<Type> >& fld = tfld();

    forAll(x1, i)
    {
        fld[i] = this->dimIntegrate(x1[i], x2[i]);
    }
    return tfld;
}


// * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * * * //

#include "DataEntryIO.C"


// ************************************************************************* //
