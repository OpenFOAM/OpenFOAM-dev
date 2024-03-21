/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "lduInterfaceFieldPtrsList.H"
#include "LduInterfaceFieldPtrsList.H"
#include "GeometricBoundaryField.H"
#include "GeometricFieldSources.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class dictionary;

// Forward declaration of friend functions and operators

template<class Type, template<class> class PatchField, class GeoMesh>
class GeometricField;

template<class Type, template<class> class PatchField, class GeoMesh>
Ostream& operator<<
(
    Ostream&,
    const GeometricField<Type, PatchField, GeoMesh>&
);

template<class Type, template<class> class PatchField, class GeoMesh>
Ostream& operator<<
(
    Ostream&,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>&
);


/*---------------------------------------------------------------------------*\
                       Class GeometricField Declaration
\*---------------------------------------------------------------------------*/

template<class Type, template<class> class PatchField, class GeoMesh>
class GeometricField
:
    public DimensionedField<Type, GeoMesh>,
    public OldTimeField<GeometricField<Type, PatchField, GeoMesh>>
{
    // Private Member Functions

        //- Read from file if it is present
        bool readIfPresent();


public:

    // Public Typedefs

        //- Type of mesh on which this GeometricField is instantiated
        typedef typename GeoMesh::Mesh Mesh;

        //- Type of the internal field from which this GeometricField is derived
        typedef DimensionedField<Type, GeoMesh> Internal;

        //- Type of the field from which this GeometricField is derived
        typedef DimensionedField<Type, GeoMesh> Base;

        //- Type of the patch field of which the Boundary is composed
        typedef PatchField<Type> Patch;

        //- Type of the boundary field
        typedef GeometricBoundaryField<Type, PatchField, GeoMesh> Boundary;

        //- Type of the field source of which the Sources is composed
        typedef typename GeoMesh::template FieldSource<Type> Source;

        //- Type of the field sources
        typedef GeometricFieldSources<Type, GeoMesh> Sources;

        //- Component type of the elements of the field
        typedef typename Field<Type>::cmptType cmptType;


private:

    // Private Data

        //-  Pointer to previous iteration (used for under-relaxation)
        mutable GeometricField<Type, PatchField, GeoMesh>* fieldPrevIterPtr_;

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


    // Static Member Functions

        //- Return a null geometric field
        inline static const GeometricField<Type, PatchField, GeoMesh>& null();


    // Constructors

        //- Constructor given IOobject, mesh, dimensions and patch field type.
        //  This allocates storage for the field but not values.
        //  Used only within this class to create TEMPORARY variables
        GeometricField
        (
            const IOobject&,
            const Mesh&,
            const dimensionSet&,
            const word& patchFieldType = PatchField<Type>::calculatedType()
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
            const word& patchFieldType = PatchField<Type>::calculatedType()
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
            const PtrList<PatchField<Type>>&,
            const HashPtrTable<Source>& = HashPtrTable<Source>()
        );

        //- Constructor from components
        GeometricField
        (
            const IOobject&,
            const Mesh&,
            const dimensionSet&,
            const Field<Type>&,
            const PtrList<PatchField<Type>>&,
            const HashPtrTable<Source>& = HashPtrTable<Source>()
        );

        //- Construct and read given IOobject
        GeometricField
        (
            const IOobject&,
            const Mesh&
        );

        //- Construct from dictionary
        GeometricField
        (
            const IOobject&,
            const Mesh&,
            const dictionary&
        );

        //- Copy constructor
        GeometricField(const GeometricField<Type, PatchField, GeoMesh>&);

        //- Move constructor
        GeometricField(GeometricField<Type, PatchField, GeoMesh>&&);

        //- Construct as copy of tmp
        GeometricField(const tmp<GeometricField<Type, PatchField, GeoMesh>>&);

        //- Construct as copy resetting IO parameters
        GeometricField
        (
            const IOobject&,
            const GeometricField<Type, PatchField, GeoMesh>&
        );

        //- Construct as copy of tmp resetting IO parameters
        GeometricField
        (
            const IOobject&,
            const tmp<GeometricField<Type, PatchField, GeoMesh>>&
        );

        //- Construct as copy resetting name
        GeometricField
        (
            const word& newName,
            const GeometricField<Type, PatchField, GeoMesh>&
        );

        //- Construct as copy of tmp resetting name
        GeometricField
        (
            const word& newName,
            const tmp<GeometricField<Type, PatchField, GeoMesh>>&
        );

        //- Construct as copy resetting IO parameters and patch field type
        GeometricField
        (
            const IOobject&,
            const GeometricField<Type, PatchField, GeoMesh>&,
            const word& patchFieldType
        );

        //- Construct as copy resetting IO parameters and boundary types
        GeometricField
        (
            const IOobject&,
            const GeometricField<Type, PatchField, GeoMesh>&,
            const wordList& patchFieldTypes,
            const wordList& actualPatchTypes = wordList(),
            const HashTable<word>& fieldSourceTypes = HashTable<word>()
        );

        //- Construct as copy of tmp resetting IO parameters and boundary types
        GeometricField
        (
            const IOobject&,
            const tmp<GeometricField<Type, PatchField, GeoMesh>>&,
            const wordList& patchFieldTypes,
            const wordList& actualPatchTypes = wordList(),
            const HashTable<word>& fieldSourceTypes = HashTable<word>()
        );

        //- Clone
        tmp<GeometricField<Type, PatchField, GeoMesh>> clone() const;

        //- Clone un-sliced
        tmp<GeometricField<Type, PatchField, GeoMesh>> cloneUnSliced() const;

        //- Return a temporary field constructed from name,
        //  internal field and list of patch fields
        static tmp<GeometricField<Type, PatchField, GeoMesh>> New
        (
            const word& name,
            const Internal&,
            const PtrList<PatchField<Type>>&,
            const HashPtrTable<Source>& = HashPtrTable<Source>()
        );

        //- Return a temporary field constructed from name, mesh, dimensionSet
        //  and patch field type.
        static tmp<GeometricField<Type, PatchField, GeoMesh>> New
        (
            const word& name,
            const Mesh&,
            const dimensionSet&,
            const word& patchFieldType = PatchField<Type>::calculatedType()
        );

        //- Return a temporary field constructed from mesh, dimensioned<Type>
        //  and patch field type.
        static tmp<GeometricField<Type, PatchField, GeoMesh>> New
        (
            const word& name,
            const Mesh&,
            const dimensioned<Type>&,
            const word& patchFieldType = PatchField<Type>::calculatedType()
        );

        //- Return a temporary field constructed from mesh, dimensioned<Type>
        //  and patch field types.
        static tmp<GeometricField<Type, PatchField, GeoMesh>> New
        (
            const word& name,
            const Mesh&,
            const dimensioned<Type>&,
            const wordList& patchFieldTypes,
            const wordList& actualPatchTypes = wordList(),
            const HashTable<word>& fieldSourceTypes = HashTable<word>()
        );

        //- Rename temporary field and return
        static tmp<GeometricField<Type, PatchField, GeoMesh>> New
        (
            const word& newName,
            const tmp<GeometricField<Type, PatchField, GeoMesh>>&
        );

        //- Rename temporary field, reset patch field type and return
        static tmp<GeometricField<Type, PatchField, GeoMesh>> New
        (
            const word& newName,
            const tmp<GeometricField<Type, PatchField, GeoMesh>>&,
            const word& patchFieldType
        );

        //- Rename and reset patch fields types of temporary field and return
        static tmp<GeometricField<Type, PatchField, GeoMesh>> New
        (
            const word& newName,
            const tmp<GeometricField<Type, PatchField, GeoMesh>>&,
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

        //- Return const-reference to the boundary field
        inline const Boundary& boundaryField() const;

        //- Return const-reference to the sources
        inline const Sources& sources() const;

        //- Store the field as the previous iteration value
        void storePrevIter() const;

        //- Delete the previous iteration field
        void clearPrevIter();

        //- Return previous iteration field
        const GeometricField<Type, PatchField, GeoMesh>& prevIter() const;

        //- Use old-time methods from the base class
        USING_OLD_TIME_FIELD(GeometricField);

        //- Correct boundary field
        void correctBoundaryConditions();

        //- Reset the field contents to the given field
        //  Used for mesh to mesh mapping
        void reset(const tmp<GeometricField<Type, PatchField, GeoMesh>>&);

        //- Does the field need a reference level for solution
        bool needReference() const;

        //- Return a component of the field
        tmp<GeometricField<cmptType, PatchField, GeoMesh>> component
        (
            const direction
        ) const;

        //- WriteData member function required by regIOobject
        bool writeData(Ostream&) const;

        //- Return transpose (only if it is a tensor field)
        tmp<GeometricField<Type, PatchField, GeoMesh>> T() const;

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
        void relax
        (
            const tmp<GeometricField<Type, PatchField, GeoMesh>>&,
            const scalar alpha
        );

        //- Relax given field with respect to the current field
        //  and reset the field to the result
        //  Relaxation factor is read from fvSolution
        void relax(const tmp<GeometricField<Type, PatchField, GeoMesh>>&);

        //- Select the final iteration parameters if `final' is true
        //  by returning the field name + "Final"
        //  otherwise the standard parameters by returning the field name
        word select(bool final) const;

        //- Helper function to write the min and max to an Ostream
        void writeMinMax(Ostream& os) const;


    // Member function *this operators

        void negate();

        void replace
        (
            const direction,
            const GeometricField<cmptType, PatchField, GeoMesh>&
        );

        void replace
        (
            const direction,
            const dimensioned<cmptType>&
        );

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

        void operator=(const GeometricField<Type, PatchField, GeoMesh>&);
        void operator=(GeometricField<Type, PatchField, GeoMesh>&&);
        void operator=(const tmp<GeometricField<Type, PatchField, GeoMesh>>&);
        void operator=(const dimensioned<Type>&);
        void operator=(const zero&);

        void operator==(const tmp<GeometricField<Type, PatchField, GeoMesh>>&);
        void operator==(const dimensioned<Type>&);
        void operator==(const zero&);

        void operator+=(const GeometricField<Type, PatchField, GeoMesh>&);
        void operator+=(const tmp<GeometricField<Type, PatchField, GeoMesh>>&);

        void operator-=(const GeometricField<Type, PatchField, GeoMesh>&);
        void operator-=(const tmp<GeometricField<Type, PatchField, GeoMesh>>&);

        void operator*=(const GeometricField<scalar, PatchField, GeoMesh>&);
        void operator*=(const tmp<GeometricField<scalar,PatchField,GeoMesh>>&);

        void operator/=(const GeometricField<scalar, PatchField, GeoMesh>&);
        void operator/=(const tmp<GeometricField<scalar,PatchField,GeoMesh>>&);

        void operator+=(const dimensioned<Type>&);
        void operator-=(const dimensioned<Type>&);

        void operator*=(const dimensioned<scalar>&);
        void operator/=(const dimensioned<scalar>&);


    // Ostream operators

        friend Ostream& operator<< <Type, PatchField, GeoMesh>
        (
            Ostream&,
            const GeometricField<Type, PatchField, GeoMesh>&
        );

        friend Ostream& operator<< <Type, PatchField, GeoMesh>
        (
            Ostream&,
            const tmp<GeometricField<Type, PatchField, GeoMesh>>&
        );
};


template<class Type, template<class> class PatchField, class GeoMesh>
Ostream& operator<<
(
    Ostream&,
    const GeometricBoundaryField<Type, PatchField, GeoMesh>&
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
