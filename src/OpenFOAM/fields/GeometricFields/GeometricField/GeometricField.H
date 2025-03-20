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
    Foam::GeometricField

Description
    Generic GeometricField class.

SourceFiles
    GeometricFieldI.H
    GeometricField.C
    GeometricFieldFunctions.H
    GeometricFieldFunctions.C

\*---------------------------------------------------------------------------*/

#ifndef GeometricField_H
#define GeometricField_H

#include "GeometricFieldFwd.H"
#include "GeometricBoundaryField.H"
#include "GeometricFieldSources.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class dictionary;

// Forward declaration of friend functions and operators

template<class Type, class GeoMesh, template<class> class PrimitiveField>
Ostream& operator<<
(
    Ostream&,
    const GeometricField<Type, GeoMesh, PrimitiveField>&
);

template<class Type, class GeoMesh, template<class> class PrimitiveField>
Ostream& operator<<
(
    Ostream&,
    const tmp<GeometricField<Type, GeoMesh, PrimitiveField>>&
);


/*---------------------------------------------------------------------------*\
                       Class GeometricField Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class GeoMesh, template<class> class PrimitiveField>
class GeometricField
:
    public DimensionedField<Type, GeoMesh, PrimitiveField>,
    public OldTimeField<GeometricField<Type, GeoMesh, PrimitiveField>>
{
    // Private Member Functions

        //- Read from file if it is present
        bool readIfPresent();


public:

    // Public Typedefs

        //- Type of mesh on which this GeometricField is instantiated
        typedef typename GeoMesh::Mesh Mesh;

        //- Type of the internal field from which this GeometricField is derived
        typedef DimensionedField<Type, GeoMesh, PrimitiveField> Internal;

        //- Type of the field from which this GeometricField is derived
        typedef DimensionedField<Type, GeoMesh, PrimitiveField> Base;

        //- Type of the patch field of which the Boundary is composed
        typedef typename GeoMesh::template PatchField<Type> Patch;

        //- Type of the boundary field
        typedef GeometricBoundaryField<Type, GeoMesh, PrimitiveField> Boundary;

        //- Type of the field source of which the Sources is composed
        typedef typename GeoMesh::template FieldSource<Type> Source;

        //- Type of the field sources
        typedef GeometricFieldSources<Type, GeoMesh, PrimitiveField> Sources;

        //- Component type of the elements of the field
        typedef typename Field<Type>::cmptType cmptType;


private:

    // Private Data

        //-  Pointer to previous iteration (used for under-relaxation)
        mutable GeometricField<Type, GeoMesh, Field>* fieldPrevIterPtr_;

        //- Boundary Type field containing boundary field values
        Boundary boundaryField_;

        //- Type field sources containing field source values
        Sources sources_;


    // Private Member Functions

        //- Read the field from the dictionary
        void readFields(const dictionary&);

        //- Read the field - create the field dictionary on-the-fly
        void readFields();


public:

    //- Runtime type information
    TypeName("GeometricField");


    //- Declare friendship with other geometric fields
    template<class Type2, class GeoMesh2, template<class> class PrimitiveField2>
    friend class GeometricField;


    // Static Member Functions

        //- Return a null geometric field
        inline static const
            GeometricField<Type, GeoMesh, PrimitiveField>& null();


    // Constructors

        //- Constructor given IOobject, mesh, dimensions and patch field type.
        //  This allocates storage for the field but not values.
        //  Used only within this class to create TEMPORARY variables
        GeometricField
        (
            const IOobject&,
            const Mesh&,
            const dimensionSet&,
            const word& patchFieldType = Patch::calculatedType()
        );

        //- Constructor given IOobject, mesh, dimensions and patch field types.
        //  This allocates storage for the field but not values.
        //  Used only within this class to create TEMPORARY variables
        GeometricField
        (
            const IOobject&,
            const Mesh&,
            const dimensionSet&,
            const wordList& wantedPatchTypes,
            const wordList& actualPatchTypes = wordList(),
            const HashTable<word>& fieldSourceTypes = HashTable<word>()
        );

        //- Constructor given IOobject, mesh, dimensioned<Type>
        //  and patch field type.
        GeometricField
        (
            const IOobject&,
            const Mesh&,
            const dimensioned<Type>&,
            const word& patchFieldType = Patch::calculatedType()
        );

        //- Constructor given IOobject, mesh, dimensioned<Type>
        //  and patch field types.
        GeometricField
        (
            const IOobject&,
            const Mesh&,
            const dimensioned<Type>&,
            const wordList& wantedPatchTypes,
            const wordList& actualPatchTypes = wordList(),
            const HashTable<word>& fieldSourceTypes = HashTable<word>()
        );

        //- Constructor from components
        GeometricField
        (
            const IOobject&,
            const Internal&,
            const PtrList<Patch>&,
            const HashPtrTable<Source>& = HashPtrTable<Source>()
        );

        //- Constructor from components
        GeometricField
        (
            const IOobject&,
            const Mesh&,
            const dimensionSet&,
            const PrimitiveField<Type>&,
            const PtrList<Patch>&,
            const HashPtrTable<Source>& = HashPtrTable<Source>()
        );

        //- Construct and read given IOobject
        GeometricField(const IOobject&, const Mesh&);

        //- Construct from dictionary
        GeometricField(const IOobject&, const Mesh&, const dictionary&);

        //- Copy constructor
        GeometricField(const GeometricField<Type, GeoMesh, PrimitiveField>&);

        //- Copy constructor
        template<template<class> class PrimitiveField2>
        explicit GeometricField
        (
            const GeometricField<Type, GeoMesh, PrimitiveField2>&
        );

        //- Move constructor
        GeometricField(GeometricField<Type, GeoMesh, PrimitiveField>&&);

        //- Construct as copy of tmp
        GeometricField
        (
            const tmp<GeometricField<Type, GeoMesh, PrimitiveField>>&
        );

        //- Construct as copy resetting IO parameters
        template<template<class> class PrimitiveField2>
        GeometricField
        (
            const IOobject&,
            const GeometricField<Type, GeoMesh, PrimitiveField2>&
        );

        //- Construct as copy of tmp resetting IO parameters
        GeometricField
        (
            const IOobject&,
            const tmp<GeometricField<Type, GeoMesh, PrimitiveField>>&
        );

        //- Construct as copy resetting name
        template<template<class> class PrimitiveField2>
        GeometricField
        (
            const word& newName,
            const GeometricField<Type, GeoMesh, PrimitiveField2>&
        );

        //- Construct as copy of tmp resetting name
        GeometricField
        (
            const word& newName,
            const tmp<GeometricField<Type, GeoMesh, PrimitiveField>>&
        );

        //- Construct as copy resetting IO parameters and patch field type
        template<template<class> class PrimitiveField2>
        GeometricField
        (
            const IOobject&,
            const GeometricField<Type, GeoMesh, PrimitiveField2>&,
            const word& patchFieldType
        );

        //- Construct as copy of tmp resetting IO parameters
        //  and patch field type
        GeometricField
        (
            const IOobject&,
            const tmp<GeometricField<Type, GeoMesh, PrimitiveField>>&,
            const word& patchFieldType
        );

        //- Constructor from internal field and patch field type
        template<template<class> class PrimitiveField2>
        GeometricField
        (
            const IOobject&,
            const DimensionedField<Type, GeoMesh, PrimitiveField2>&,
            const word& patchFieldType
        );

        //- Constructor from temporary internal field and patch field type
        GeometricField
        (
            const IOobject&,
            const tmp<Internal>&,
            const word& patchFieldType
        );

        //- Construct as copy resetting IO parameters and patch field types
        template<template<class> class PrimitiveField2>
        GeometricField
        (
            const IOobject&,
            const GeometricField<Type, GeoMesh, PrimitiveField2>&,
            const wordList& patchFieldTypes,
            const wordList& actualPatchTypes = wordList(),
            const HashTable<word>& fieldSourceTypes = HashTable<word>()
        );

        //- Construct as copy of tmp resetting IO parameters
        //  and patch field types
        GeometricField
        (
            const IOobject&,
            const tmp<GeometricField<Type, GeoMesh, PrimitiveField>>&,
            const wordList& patchFieldTypes,
            const wordList& actualPatchTypes = wordList(),
            const HashTable<word>& fieldSourceTypes = HashTable<word>()
        );

        //- Construct from IOobject, internal field and patch field types
        template<template<class> class PrimitiveField2>
        GeometricField
        (
            const IOobject&,
            const DimensionedField<Type, GeoMesh, PrimitiveField2>&,
            const wordList& patchFieldTypes,
            const wordList& actualPatchTypes = wordList(),
            const HashTable<word>& fieldSourceTypes = HashTable<word>()
        );

        //- Construct from IOobject, tmp internal field and patch field types
        GeometricField
        (
            const IOobject&,
            const tmp<Internal>&,
            const wordList& patchFieldTypes,
            const wordList& actualPatchTypes = wordList(),
            const HashTable<word>& fieldSourceTypes = HashTable<word>()
        );

        //- Clone
        tmp<GeometricField<Type, GeoMesh, PrimitiveField>> clone() const;

        //- Clone un-sliced
        tmp<GeometricField<Type, GeoMesh, PrimitiveField>>
            cloneUnSliced() const;

        //- Return a temporary field constructed from name,
        //  internal field and list of patch fields
        static tmp<GeometricField<Type, GeoMesh, PrimitiveField>> New
        (
            const word& name,
            const Internal&,
            const PtrList<Patch>&,
            const HashPtrTable<Source>& = HashPtrTable<Source>()
        );

        //- Return a temporary field constructed from name, mesh, dimensionSet
        //  and patch field type
        static tmp<GeometricField<Type, GeoMesh, PrimitiveField>> New
        (
            const word& name,
            const Mesh&,
            const dimensionSet&,
            const word& patchFieldType = Patch::calculatedType()
        );

        //- Return a temporary field constructed from mesh, dimensioned<Type>
        //  and patch field type
        static tmp<GeometricField<Type, GeoMesh, PrimitiveField>> New
        (
            const word& name,
            const Mesh&,
            const dimensioned<Type>&,
            const word& patchFieldType = Patch::calculatedType()
        );

        //- Return a temporary field constructed from mesh, dimensioned<Type>
        //  and patch field types
        static tmp<GeometricField<Type, GeoMesh, PrimitiveField>> New
        (
            const word& name,
            const Mesh&,
            const dimensioned<Type>&,
            const wordList& patchFieldTypes,
            const wordList& actualPatchTypes = wordList(),
            const HashTable<word>& fieldSourceTypes = HashTable<word>()
        );

        //- Rename temporary field and return
        static tmp<GeometricField<Type, GeoMesh, PrimitiveField>> New
        (
            const word& newName,
            const tmp<GeometricField<Type, GeoMesh, PrimitiveField>>&
        );

        //- Rename field, reset patch field type and return
        template<template<class> class PrimitiveField2>
        static tmp<GeometricField<Type, GeoMesh, PrimitiveField>> New
        (
            const word& newName,
            const GeometricField<Type, GeoMesh, PrimitiveField2>&,
            const word& patchFieldType
        );

        //- Rename temporary field, reset patch field type and return
        static tmp<GeometricField<Type, GeoMesh, PrimitiveField>> New
        (
            const word& newName,
            const tmp<GeometricField<Type, GeoMesh, PrimitiveField>>&,
            const word& patchFieldType
        );

        //- Return a temporary field constructed from name,
        //  internal field and patch field type
        template<template<class> class PrimitiveField2>
        static tmp<GeometricField<Type, GeoMesh, PrimitiveField>> New
        (
            const word& name,
            const DimensionedField<Type, GeoMesh, PrimitiveField2>&,
            const word& patchFieldType
        );

        //- Return a temporary field constructed from name,
        //  temporary internal field and patch field type
        static tmp<GeometricField<Type, GeoMesh, PrimitiveField>> New
        (
            const word& name,
            const tmp<Internal>&,
            const word& patchFieldType
        );

        //- Rename and reset patch fields types of field and return
        template<template<class> class PrimitiveField2>
        static tmp<GeometricField<Type, GeoMesh, PrimitiveField>> New
        (
            const word& newName,
            const GeometricField<Type, GeoMesh, PrimitiveField2>&,
            const wordList& patchFieldTypes,
            const wordList& actualPatchTypes = wordList(),
            const HashTable<word>& fieldSourceTypes = HashTable<word>()
        );

        //- Rename and reset patch fields types of temporary field and return
        static tmp<GeometricField<Type, GeoMesh, PrimitiveField>> New
        (
            const word& newName,
            const tmp<GeometricField<Type, GeoMesh, PrimitiveField>>&,
            const wordList& patchFieldTypes,
            const wordList& actualPatchTypes = wordList(),
            const HashTable<word>& fieldSourceTypes = HashTable<word>()
        );

        //- Return a temporary field constructed from IOobject,
        //  internal field and patch field types
        template<template<class> class PrimitiveField2>
        static tmp<GeometricField<Type, GeoMesh, PrimitiveField>> New
        (
            const word& newName,
            const DimensionedField<Type, GeoMesh, PrimitiveField2>&,
            const wordList& patchFieldTypes,
            const wordList& actualPatchTypes = wordList(),
            const HashTable<word>& fieldSourceTypes = HashTable<word>()
        );

        //- Return a temporary field constructed from IOobject,
        //  temporary internal field and patch field types
        static tmp<GeometricField<Type, GeoMesh, PrimitiveField>> New
        (
            const word& newName,
            const tmp<Internal>&,
            const wordList& patchFieldTypes,
            const wordList& actualPatchTypes = wordList(),
            const HashTable<word>& fieldSourceTypes = HashTable<word>()
        );


    //- Destructor
    virtual ~GeometricField();


    // Member Functions

        //- Return a reference to the dimensioned internal field
        //  Note: this increments the event counter and checks the
        //  old-time fields; avoid in loops.
        Internal& internalFieldRef();

        //- Return a const-reference to the dimensioned internal field
        inline const Internal& internalField() const;

        //- Return a const-reference to the dimensioned internal field
        //  of a "vol" field.  Useful in the formulation of source-terms
        //  for FV equations
        inline const Internal& v() const;

        //- Return a reference to the primitive field
        //  Note: this increments the event counter and checks the
        //  old-time fields; avoid in loops.
        typename Internal::FieldType& primitiveFieldRef();

        //- Return a const-reference to the primitive field
        inline const typename Internal::FieldType& primitiveField() const;

        //- Return a reference to the boundary field
        //  Note: this increments the event counter and checks the
        //  old-time fields; avoid in loops.
        Boundary& boundaryFieldRef();

        //- Return a reference to the boundary field without storing old times
        //  Note: this increments the event counter; avoid in loops.
        Boundary& boundaryFieldRefNoStoreOldTimes();

        //- Return const-reference to the boundary field
        inline const Boundary& boundaryField() const;

        //- Return a reference to the sources
        //  Note: this increments the event counter and checks the
        //  old-time fields; avoid in loops.
        Sources& sourcesRef();

        //- Return const-reference to the sources
        inline const Sources& sources() const;

        //- Store the field as the previous iteration value
        void storePrevIter() const;

        //- Delete the previous iteration field
        void clearPrevIter();

        //- Return previous iteration field
        const GeometricField<Type, GeoMesh, Field>& prevIter() const;

        //- Use old-time methods from the base class
        USING_OLD_TIME_FIELD(GeometricField);

        //- Correct boundary field
        void correctBoundaryConditions();

        //- Reset the field contents to the given field
        //  Used for mesh to mesh mapping
        template<template<class> class PrimitiveField2>
        void reset(const GeometricField<Type, GeoMesh, PrimitiveField2>&);

        //- Reset the field contents to the given tmp-field
        //  Used for mesh to mesh mapping
        void reset(const tmp<GeometricField<Type, GeoMesh, PrimitiveField>>&);

        //- Reset the field contents to the given tmp-field
        //  Used for mesh to mesh mapping
        template<template<class> class PrimitiveField2>
        void reset(const tmp<GeometricField<Type, GeoMesh, PrimitiveField2>>&);

        //- Does the field need a reference level for solution
        bool needReference() const;

        //- Return a component of the field
        tmp<GeometricField<cmptType, GeoMesh, Field>> component
        (
            const direction
        ) const;

        //- WriteData member function required by regIOobject
        bool writeData(Ostream&) const;

        //- Return transpose (only if it is a tensor field)
        tmp<GeometricField<Type, GeoMesh, Field>> T() const;

        //- Relax field (for steady-state solution).
        //  alpha >= 1 : no relaxation
        //  alpha < 1  : relaxation
        void relax(const scalar alpha);

        //- Return the field relaxation factor read from fvSolution
        //  or 1 if not specified
        scalar relaxationFactor() const;

        //- Relax current field with respect to the cached previous iteration.
        //  Relaxation factor is read from fvSolution
        void relax();

        //- Relax given field with respect to the current field
        //  and reset the field to the result
        template<template<class> class PrimitiveField2>
        void relax
        (
            const tmp<GeometricField<Type, GeoMesh, PrimitiveField2>>&,
            const scalar alpha
        );

        //- Relax given field with respect to the current field
        //  and reset the field to the result
        //  Relaxation factor is read from fvSolution
        template<template<class> class PrimitiveField2>
        void relax(const tmp<GeometricField<Type, GeoMesh, PrimitiveField2>>&);

        //- Select the final iteration parameters if `final' is true
        //  by returning the field name + "Final"
        //  otherwise the standard parameters by returning the field name
        word select(bool final) const;

        //- Helper function to write the min and max to an Ostream
        void writeMinMax(Ostream& os) const;


    // Member function *this operators

        void negate();

        //- Replace a component field of the field
        template<template<class> class PrimitiveField2>
        void replace
        (
            const direction,
            const GeometricField<cmptType, GeoMesh, PrimitiveField2>&
        );

        //- Replace a component field of the field
        template<template<class> class PrimitiveField2>
        void replace
        (
            const direction,
            const tmp<GeometricField<cmptType, GeoMesh, PrimitiveField2>>&
        );

        //- Replace a component field of the field
        void replace(const direction, const dimensioned<cmptType>&);

        void max(const dimensioned<Type>&);

        void min(const dimensioned<Type>&);

        void maxMin
        (
            const dimensioned<Type>& minDt,
            const dimensioned<Type>& maxDt
        );


    // Member Operators

        //- Return a const-reference to the dimensioned internal field
        //  Useful in the formulation of source-terms for FV equations
        inline const Internal& operator()() const;

        void operator=
        (
            const GeometricField<Type, GeoMesh, PrimitiveField>&
        );
        void operator=
        (
            GeometricField<Type, GeoMesh, PrimitiveField>&&
        );
        template<template<class> class PrimitiveField2>
        void operator=
        (
            const GeometricField<Type, GeoMesh, PrimitiveField2>&
        );
        void operator=
        (
            const tmp<GeometricField<Type, GeoMesh, PrimitiveField>>&
        );
        template<template<class> class PrimitiveField2>
        void operator=
        (
            const tmp<GeometricField<Type, GeoMesh, PrimitiveField2>>&
        );
        void operator=(const dimensioned<Type>&);
        void operator=(const zero&);

        template<template<class> class PrimitiveField2>
        void operator==
        (
            const GeometricField<Type, GeoMesh, PrimitiveField2>&
        );
        void operator==
        (
            const tmp<GeometricField<Type, GeoMesh, PrimitiveField>>&
        );
        template<template<class> class PrimitiveField2>
        void operator==
        (
            const tmp<GeometricField<Type, GeoMesh, PrimitiveField2>>&
        );
        void operator==(const dimensioned<Type>&);
        void operator==(const zero&);

        template<template<class> class PrimitiveField2>
        void operator+=
        (
            const GeometricField<Type, GeoMesh, PrimitiveField2>&
        );
        template<template<class> class PrimitiveField2>
        void operator+=
        (
            const tmp<GeometricField<Type, GeoMesh, PrimitiveField2>>&
        );
        void operator+=(const dimensioned<Type>&);

        template<template<class> class PrimitiveField2>
        void operator-=
        (
            const GeometricField<Type, GeoMesh, PrimitiveField2>&
        );
        template<template<class> class PrimitiveField2>
        void operator-=
        (
            const tmp<GeometricField<Type, GeoMesh, PrimitiveField2>>&
        );
        void operator-=(const dimensioned<Type>&);

        template<template<class> class PrimitiveField2>
        void operator*=
        (
            const GeometricField<scalar, GeoMesh, PrimitiveField2>&
        );
        template<template<class> class PrimitiveField2>
        void operator*=
        (
            const tmp<GeometricField<scalar, GeoMesh, PrimitiveField2>>&
        );
        void operator*=(const dimensioned<scalar>&);

        template<template<class> class PrimitiveField2>
        void operator/=
        (
            const GeometricField<scalar, GeoMesh, PrimitiveField2>&
        );
        template<template<class> class PrimitiveField2>
        void operator/=
        (
            const tmp<GeometricField<scalar, GeoMesh, PrimitiveField2>>&
        );
        void operator/=(const dimensioned<scalar>&);


    // Ostream operators

        friend Ostream& operator<< <Type, GeoMesh, PrimitiveField>
        (
            Ostream&,
            const GeometricField<Type, GeoMesh, PrimitiveField>&
        );

        friend Ostream& operator<< <Type, GeoMesh, PrimitiveField>
        (
            Ostream&,
            const tmp<GeometricField<Type, GeoMesh, PrimitiveField>>&
        );
};


template<class Type, class GeoMesh, template<class> class PrimitiveField>
Ostream& operator<<
(
    Ostream&,
    const GeometricBoundaryField<Type, GeoMesh, PrimitiveField>&
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "GeometricFieldI.H"

#ifdef NoRepository
    #include "GeometricField.C"
#endif

#include "GeometricFieldFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
