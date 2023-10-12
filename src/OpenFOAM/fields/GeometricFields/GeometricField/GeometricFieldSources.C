/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "GeometricFieldSources.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
Foam::GeometricFieldSources<Type, GeoMesh>::GeometricFieldSources()
:
    HashPtrTable<Source>(),
    errorLocation_()
{}


template<class Type, class GeoMesh>
Foam::GeometricFieldSources<Type, GeoMesh>::GeometricFieldSources
(
    const DimensionedField<Type, GeoMesh>& iF,
    const HashPtrTable<Source>& mtf
)
:
    HashPtrTable<Source>(mtf.capacity()),
    errorLocation_()
{
    forAllConstIter(typename HashPtrTable<Source>, mtf, iter)
    {
        this->set(iter.key(), iter()->clone(iF).ptr());
    }
}


template<class Type, class GeoMesh>
Foam::GeometricFieldSources<Type, GeoMesh>::GeometricFieldSources
(
    const DimensionedField<Type, GeoMesh>& iF,
    const GeometricFieldSources<Type, GeoMesh>& mtf
)
:
    GeometricFieldSources(iF, static_cast<const HashPtrTable<Source>&>(mtf))
{}


template<class Type, class GeoMesh>
Foam::GeometricFieldSources<Type, GeoMesh>::GeometricFieldSources
(
    const DimensionedField<Type, GeoMesh>& iF,
    const dictionary& dict
)
:
    HashPtrTable<Source>(),
    errorLocation_()
{
    readField(iF, dict);
}


template<class Type, class GeoMesh>
Foam::GeometricFieldSources<Type, GeoMesh>::GeometricFieldSources
(
    const DimensionedField<Type, GeoMesh>& iF,
    const HashTable<word>& types
)
:
    HashPtrTable<Source>(),
    errorLocation_()
{
    forAllConstIter(typename HashTable<word>, types, iter)
    {
        this->set
        (
            iter.key(),
            Source::New(iter(), iF).ptr()
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
Foam::GeometricFieldSources<Type, GeoMesh>::
~GeometricFieldSources()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
const Foam::HashPtrTable
<
    typename Foam::GeometricFieldSources<Type, GeoMesh>::Source
>&
Foam::GeometricFieldSources<Type, GeoMesh>::table() const
{
    return *this;
}


template<class Type, class GeoMesh>
Foam::HashPtrTable
<
    typename Foam::GeometricFieldSources<Type, GeoMesh>::Source
>&
Foam::GeometricFieldSources<Type, GeoMesh>::table()
{
    return *this;
}


template<class Type, class GeoMesh>
Foam::HashTable<Foam::word>
Foam::GeometricFieldSources<Type, GeoMesh>::types() const
{
    HashTable<word> result;

    forAllConstIter(typename HashPtrTable<Source>, *this, iter)
    {
        result.insert(iter.key(), iter()->type());
    }

    return result;
}


template<class Type, class GeoMesh>
void Foam::GeometricFieldSources<Type, GeoMesh>::readField
(
    const DimensionedField<Type, GeoMesh>& field,
    const dictionary& dict
)
{
    this->clear();

    errorLocation_ = IOerrorLocation(dict);

    forAllConstIter(dictionary, dict, iter)
    {
        if (iter().isDict())
        {
            this->set
            (
                iter().keyword(),
                Source::New(field, iter().dict()).ptr()
            );
        }
    }
}


template<class Type, class GeoMesh>
void Foam::GeometricFieldSources<Type, GeoMesh>::reset
(
    const GeometricFieldSources<Type, GeoMesh>& mtf
)
{
    this->clear();

    if (mtf.empty()) return;

    const DimensionedField<Type, GeoMesh>& iF =
        (**mtf.HashPtrTable<Source>::begin()).internalField();

    errorLocation_ = mtf.errorLocation_;

    forAllConstIter(typename HashPtrTable<Source>, mtf, iter)
    {
        this->set(iter.key(), iter()->clone(iF).ptr());
    }
}


template<class Type, class GeoMesh>
void Foam::GeometricFieldSources<Type, GeoMesh>::writeEntry
(
    const word& keyword,
    Ostream& os
) const
{
    os  << keyword << nl << token::BEGIN_BLOCK << incrIndent << nl;

    forAllConstIter(typename HashPtrTable<Source>, *this, iter)
    {
        os  << indent << iter.key() << nl
            << indent << token::BEGIN_BLOCK << nl
            << incrIndent << *iter() << decrIndent
            << indent << token::END_BLOCK << endl;
    }

    os  << decrIndent << token::END_BLOCK << endl;

    // Check state of IOstream
    os.check
    (
        "GeometricFieldSources<Type, GeoMesh>::"
        "writeEntry(const word& keyword, Ostream& os) const"
    );
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
const typename Foam::GeometricFieldSources<Type, GeoMesh>::Source&
Foam::GeometricFieldSources<Type, GeoMesh>::operator[]
(
    const word& sourceName
) const
{
    typename HashPtrTable<Source>::const_iterator iter = this->find(sourceName);

    if (iter == this->end())
    {
        FatalIOErrorInFunction(errorLocation_)
            << "Cannot find fieldSource entry for " << sourceName
            << exit(FatalIOError);
    }

    return **iter;
}


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const GeometricFieldSources<Type, GeoMesh>& bf
)
{
    typedef typename GeoMesh::template FieldSource<Type> Source;
    os << static_cast<const HashPtrTable<Source>&>(bf);
    return os;
}


// ************************************************************************* //
