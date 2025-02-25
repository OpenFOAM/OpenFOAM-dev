/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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
    Foam::DimensionedField

Description
    Field with dimensions and associated with geometry type GeoMesh which is
    used to size the field and a reference to it is maintained.

SourceFiles
    DimensionedFieldI.H
    DimensionedField.C
    DimensionedFieldIO.C

\*---------------------------------------------------------------------------*/

#ifndef DimensionedField_H
#define DimensionedField_H

#include "DimensionedFieldFwd.H"
#include "regIOobject.H"
#include "Field.H"
#include "OldTimeField.H"
#include "dimensionedType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Type, class GeoMesh, template<class> class PrimitiveField>
class DimensionedField;

template<class Type, class GeoMesh, template<class> class PrimitiveField>
Ostream& operator<<
(
    Ostream&,
    const DimensionedField<Type, GeoMesh, PrimitiveField>&
);

template<class Type, class GeoMesh, template<class> class PrimitiveField>
Ostream& operator<<
(
    Ostream&,
    const tmp<DimensionedField<Type, GeoMesh, PrimitiveField>>&
);


/*---------------------------------------------------------------------------*\
                      Class DimensionedField Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class GeoMesh, template<class> class PrimitiveField>
class DimensionedField
:
    public regIOobject,
    public PrimitiveField<Type>,
    public OldTimeField<DimensionedField<Type, GeoMesh, PrimitiveField>>
{

public:

    // Public Typedefs

        //- Type of mesh on which this DimensionedField is instantiated
        typedef typename GeoMesh::Mesh Mesh;

        //- Type of the field from which this DimensionedField is derived
        typedef PrimitiveField<Type> FieldType;

        //- Type of the field from which this DimensionedField is derived
        typedef PrimitiveField<Type> Base;

        //- Component type of the elements of the field
        typedef typename PrimitiveField<Type>::cmptType cmptType;


private:

    // Private Data

        //- Reference to mesh
        const Mesh& mesh_;

        //- Dimension set for this field
        dimensionSet dimensions_;


    // Private Member Functions

        //- Read from file if it is present
        bool readIfPresent(const word& fieldDictEntry = "value");


public:

    //- Runtime type information
    TypeName("DimensionedField");


    //- Declare friendship with other dimensioned fields
    template<class Type2, class GeoMesh2, template<class> class PrimitiveField2>
    friend class DimensionedField;


    // Static Member Functions

        //- Return a null DimensionedField
        inline static const
            DimensionedField<Type, GeoMesh, PrimitiveField>& null();


    // Constructors

        //- Construct from components
        DimensionedField
        (
            const IOobject&,
            const Mesh& mesh,
            const dimensionSet&,
            const PrimitiveField<Type>&
        );

        //- Construct from components
        //  Used for temporary fields which are initialised after construction
        DimensionedField
        (
            const IOobject&,
            const Mesh& mesh,
            const dimensionSet&,
            const bool checkIOFlags = true
        );

        //- Construct from components
        DimensionedField
        (
            const IOobject&,
            const Mesh& mesh,
            const dimensioned<Type>&,
            const bool checkIOFlags = true
        );

        //- Construct from Istream
        DimensionedField
        (
            const IOobject&,
            const Mesh& mesh,
            const word& fieldDictEntry="value"
        );

        //- Construct from dictionary
        DimensionedField
        (
            const IOobject&,
            const Mesh& mesh,
            const dictionary& fieldDict,
            const word& fieldDictEntry="value"
        );

        //- Copy constructor
        DimensionedField
        (
            const DimensionedField<Type, GeoMesh, PrimitiveField>&
        );

        //- Move constructor
        DimensionedField
        (
            DimensionedField<Type, GeoMesh, PrimitiveField>&&
        );

        //- Copy constructor
        template<template<class> class PrimitiveField2>
        explicit DimensionedField
        (
            const DimensionedField<Type, GeoMesh, PrimitiveField2>&
        );

        //- Copy constructor or reuse as specified.
        template<template<class> class PrimitiveField2>
        DimensionedField
        (
            DimensionedField<Type, GeoMesh, PrimitiveField2>&,
            const bool reuse
        );

        //- Copy constructor of tmp<DimensionedField> deleting argument
        DimensionedField
        (
            const tmp<DimensionedField<Type, GeoMesh, PrimitiveField>>&
        );

        //- Copy constructor resetting IO parameters
        template<template<class> class PrimitiveField2>
        DimensionedField
        (
            const IOobject&,
            const DimensionedField<Type, GeoMesh, PrimitiveField2>&,
            const bool checkIOFlags = true
        );

        //- Copy constructor resetting IO parameters and reuse as specified.
        template<template<class> class PrimitiveField2>
        DimensionedField
        (
            const IOobject&,
            DimensionedField<Type, GeoMesh, PrimitiveField2>&,
            const bool reuse,
            const bool checkIOFlags = true
        );

        //- Copy constructor from tmp resetting IO parameters
        DimensionedField
        (
            const IOobject&,
            const tmp<DimensionedField<Type, GeoMesh, PrimitiveField>>&,
            const bool checkIOFlags = true
        );

        //- Copy constructor resetting name
        template<template<class> class PrimitiveField2>
        DimensionedField
        (
            const word& newName,
            const DimensionedField<Type, GeoMesh, PrimitiveField2>&
        );

        //- Copy constructor resetting name and reuse as specified.
        template<template<class> class PrimitiveField2>
        DimensionedField
        (
            const word& newName,
            DimensionedField<Type, GeoMesh, PrimitiveField2>&,
            const bool reuse
        );

        //- Copy constructor from tmp resetting name
        DimensionedField
        (
            const word& newName,
            const tmp<DimensionedField<Type, GeoMesh, PrimitiveField>>&
        );

        //- Clone
        tmp<DimensionedField<Type, GeoMesh, PrimitiveField>> clone() const;

        //- Return a temporary field constructed from name, mesh,
        //  dimensionSet and field
        static tmp<DimensionedField<Type, GeoMesh, PrimitiveField>> New
        (
            const word& name,
            const Mesh& mesh,
            const dimensionSet&,
            const PrimitiveField<Type>&
        );

        //- Return a temporary field constructed from name, mesh
        //  and dimensionSet
        static tmp<DimensionedField<Type, GeoMesh, PrimitiveField>> New
        (
            const word& name,
            const Mesh& mesh,
            const dimensionSet&
        );

        //- Return a temporary field constructed from name, mesh
        //  and dimensionType providing dimensions and value
        static tmp<DimensionedField<Type, GeoMesh, PrimitiveField>> New
        (
            const word& name,
            const Mesh& mesh,
            const dimensioned<Type>&
        );

        //- Return a temporary field constructed from name and a field
        template<template<class> class PrimitiveField2>
        static tmp<DimensionedField<Type, GeoMesh, PrimitiveField>> New
        (
            const word& newName,
            const DimensionedField<Type, GeoMesh, PrimitiveField2>&
        );

        //- Return renamed temporary field
        static tmp<DimensionedField<Type, GeoMesh, PrimitiveField>> New
        (
            const word& newName,
            const tmp<DimensionedField<Type, GeoMesh, PrimitiveField>>&
        );


    //- Destructor
    virtual ~DimensionedField();


    // Member Functions

        void readField
        (
            const dictionary& fieldDict,
            const word& fieldDictEntry = "value"
        );

        //- Return mesh
        inline const Mesh& mesh() const;

        //- Return dimensions
        inline const dimensionSet& dimensions() const;

        //- Return non-const access to dimensions
        inline dimensionSet& dimensions();

        //- Return a reference to the internal field
        PrimitiveField<Type>& primitiveFieldRef();

        //- Return a const-reference to the primitive field
        inline const PrimitiveField<Type>& primitiveField() const;

        //- Use old-time methods from the base class
        USING_OLD_TIME_FIELD(DimensionedField);

        //- Return a component field of the field
        tmp<DimensionedField<cmptType, GeoMesh, Field>> component
        (
            const direction
        ) const;

        //- Replace a component field of the field
        template<template<class> class PrimitiveField2>
        void replace
        (
            const direction,
            const DimensionedField<cmptType, GeoMesh, PrimitiveField2>&
        );

        //- Replace a component field of the field
        template<template<class> class PrimitiveField2>
        void replace
        (
            const direction,
            const tmp<DimensionedField<cmptType, GeoMesh, PrimitiveField2>>&
        );

        //- Return the field transpose (only defined for second rank tensors)
        tmp<DimensionedField<Type, GeoMesh, Field>> T() const;

        //- Calculate and return arithmetic average
        dimensioned<Type> average() const;

        //- Calculate and return weighted average
        template<template<class> class PrimitiveField2>
        dimensioned<Type> weightedAverage
        (
            const DimensionedField<scalar, GeoMesh, PrimitiveField2>&
        ) const;

        //- Calculate and return weighted average
        template<template<class> class PrimitiveField2>
        dimensioned<Type> weightedAverage
        (
            const tmp<DimensionedField<scalar, GeoMesh, PrimitiveField2>>&
        ) const;

        //- Reset the field values to the given field
        //  Used for mesh to mesh mapping
        template<template<class> class PrimitiveField2>
        void reset(const DimensionedField<Type, GeoMesh, PrimitiveField2>&);

        //- Reset the field values to the given tmp-field
        //  Used for mesh to mesh mapping
        void reset(const tmp<DimensionedField<Type, GeoMesh, PrimitiveField>>&);

        //- Reset the field values to the given tmp-field
        //  Used for mesh to mesh mapping
        template<template<class> class PrimitiveField2>
        void reset
        (
            const tmp<DimensionedField<Type, GeoMesh, PrimitiveField2>>&
        );


        // Write

            bool writeData(Ostream&, const word& fieldDictEntry) const;

            bool writeData(Ostream&) const;


    // Member Operators

        void operator=(const DimensionedField<Type, GeoMesh, PrimitiveField>&);
        void operator=(DimensionedField<Type, GeoMesh, PrimitiveField>&&);
        void operator=
        (
            const tmp<DimensionedField<Type, GeoMesh, PrimitiveField>>&
        );
        template<template<class> class PrimitiveField2>
        void operator=
        (
            const DimensionedField<Type, GeoMesh, PrimitiveField2>&
        );
        template<template<class> class PrimitiveField2>
        void operator=
        (
            const tmp<DimensionedField<Type, GeoMesh, PrimitiveField2>>&
        );
        void operator=(const dimensioned<Type>&);
        void operator=(const zero&);

        template<template<class> class PrimitiveField2>
        void operator==
        (
            const DimensionedField<Type, GeoMesh, PrimitiveField2>&
        );
        template<template<class> class PrimitiveField2>
        void operator==
        (
            const tmp<DimensionedField<Type, GeoMesh, PrimitiveField2>>&
        );
        void operator==(const dimensioned<Type>&);
        void operator==(const zero&);

        template<template<class> class PrimitiveField2>
        void operator+=
        (
            const DimensionedField<Type, GeoMesh, PrimitiveField2>&
        );
        template<template<class> class PrimitiveField2>
        void operator+=
        (
            const tmp<DimensionedField<Type, GeoMesh, PrimitiveField2>>&
        );
        void operator+=(const dimensioned<Type>&);

        template<template<class> class PrimitiveField2>
        void operator-=
        (
            const DimensionedField<Type, GeoMesh, PrimitiveField2>&
        );
        template<template<class> class PrimitiveField2>
        void operator-=
        (
            const tmp<DimensionedField<Type, GeoMesh, PrimitiveField2>>&
        );
        void operator-=(const dimensioned<Type>&);

        template<template<class> class PrimitiveField2>
        void operator*=
        (
            const DimensionedField<scalar, GeoMesh, PrimitiveField2>&
        );
        template<template<class> class PrimitiveField2>
        void operator*=
        (
            const tmp<DimensionedField<scalar, GeoMesh, PrimitiveField2>>&
        );
        void operator*=(const dimensioned<scalar>&);

        template<template<class> class PrimitiveField2>
        void operator/=
        (
            const DimensionedField<scalar, GeoMesh, PrimitiveField2>&
        );
        template<template<class> class PrimitiveField2>
        void operator/=
        (
            const tmp<DimensionedField<scalar, GeoMesh, PrimitiveField2>>&
        );
        void operator/=(const dimensioned<scalar>&);


    // Ostream Operators

        friend Ostream& operator<< <Type, GeoMesh, PrimitiveField>
        (
            Ostream&,
            const DimensionedField<Type, GeoMesh, PrimitiveField>&
        );

        friend Ostream& operator<< <Type, GeoMesh, PrimitiveField>
        (
            Ostream&,
            const tmp<DimensionedField<Type, GeoMesh, PrimitiveField>>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DimensionedFieldI.H"
#include "DimensionedFieldFunctions.H"

#ifdef NoRepository
    #include "DimensionedField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
