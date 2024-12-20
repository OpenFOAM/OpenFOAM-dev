/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2024 OpenFOAM Foundation
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

Class
    Foam::GeometricFieldSources

Description
    Part of a geometric field used for setting the values associated with
    optional sources

SourceFiles
    GeometricFieldSources.C

\*---------------------------------------------------------------------------*/

#ifndef GeometricFieldSources_H
#define GeometricFieldSources_H

#include "dimensionedTypes.H"
#include "DimensionedField.H"
#include "HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class dictionary;

/*---------------------------------------------------------------------------*\
                      Class GeometricFieldSources Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class GeoMesh, template<class> class PrimitiveField>
class GeometricFieldSources
:
    private HashPtrTable<typename GeoMesh::template FieldSource<Type>>
{
public:

    // Public Typedefs

        //- Type of the internal field from which this GeometricField is derived
        typedef DimensionedField<Type, GeoMesh, PrimitiveField> Internal;

        //- Type of the field source of which this field sources is composed
        typedef typename GeoMesh::template FieldSource<Type> Source;


private:

    // Private Data

        //- Cached IO error location for delayed error messages
        IOerrorLocation errorLocation_;


public:

    //- Declare friendship with other geometric field sources
    template<class Type2, class GeoMesh2, template<class> class PrimitiveField2>
    friend class GeometricFieldSources;


    // Constructors

        //- Construct null
        GeometricFieldSources();

        //- Construct from a table of field sources
        GeometricFieldSources(const Internal&, const HashPtrTable<Source>&);

        //- Construct as copy setting the reference to the internal field
        GeometricFieldSources(const Internal&, const GeometricFieldSources&);

        //- Construct as copy setting the reference to the internal field
        template<template<class> class PrimitiveField2>
        explicit GeometricFieldSources
        (
            const Internal&,
            const GeometricFieldSources<Type, GeoMesh, PrimitiveField2>&
        );

        //- Construct from dictionary
        GeometricFieldSources(const Internal&, const dictionary&);

        //- Construct from types
        GeometricFieldSources(const Internal&, const HashTable<word>&);

        //- Copy constructor deleted
        //  as it would not set the internalField reference correctly
        GeometricFieldSources(const GeometricFieldSources&) = delete;

        //- Move constructor deleted
        //  as it would not set the internalField reference correctly
        GeometricFieldSources(GeometricFieldSources&&) = delete;


    //- Destructor
    ~GeometricFieldSources();


    // Member Functions

        //- Access the underlying field table
        const HashPtrTable<Source>& table() const;

        //- Access the underlying field table
        HashPtrTable<Source>& table();

        //- Inherit empty test
        using HashPtrTable<Source>::empty;

        //- Return a map from the source name to the field source type
        HashTable<word> types() const;

        //- Read the sources
        void readField(const Internal& field, const dictionary& dict);

        //- Reset the boundary field contents to the given field
        //  Used for mesh to mesh mapping
        void reset(const Internal&, const GeometricFieldSources&);

        //- Write sources as dictionary entry
        void writeEntry(const word& keyword, Ostream& os) const;


    // Member operators

        //- Find and return a field source
        const Source& operator[](const word& sourceName) const;
};


template<class Type, class GeoMesh, template<class> class PrimitiveField>
Ostream& operator<<
(
    Ostream&,
    const GeometricFieldSources<Type, GeoMesh, PrimitiveField>&
);


/*---------------------------------------------------------------------------*\
                      Class NoFieldSource Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class NoFieldSource
{
public:

    // Constructors

        //- Construct null
        NoFieldSource()
        {}

        //- Construct from internal field
        template<class InternalField>
        NoFieldSource(const InternalField&)
        {}

        //- Dummy clone
        template<class InternalField>
        autoPtr<NoFieldSource<Type>> clone(const InternalField&) const
        {
            return autoPtr<NoFieldSource<Type>>(new NoFieldSource<Type>());
        }


    // Selectors

        //- Dummy selector
        template<class InternalField>
        static autoPtr<NoFieldSource<Type>> New
        (
            const word&,
            const InternalField&
        )
        {
            return autoPtr<NoFieldSource<Type>>(new NoFieldSource<Type>());
        }

        //- Dummy selector
        template<class InternalField>
        static autoPtr<NoFieldSource<Type>> New
        (
            const InternalField&,
            const dictionary&
        )
        {
            return autoPtr<NoFieldSource<Type>>(new NoFieldSource<Type>());
        }


    // Member Functions

        //- Dummy return of the internal field reference
        const nil internalField() const
        {
            return nil();
        }
};


template<class Type>
Ostream& operator<<(Ostream& os, const NoFieldSource<Type>&)
{
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "GeometricFieldSources.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
