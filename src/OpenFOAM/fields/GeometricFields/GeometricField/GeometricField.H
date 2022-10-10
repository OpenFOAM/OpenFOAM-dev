/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    public DimensionedField<Type, GeoMesh>
{
    // Private Member Functions

        //- Read from file if it is present
        bool readIfPresent();

        //- Read old time field from file if it is present
        bool readOldTimeIfPresent();


public:

    // Public Typedefs

        //- Type of mesh on which this GeometricField is instantiated
        typedef typename GeoMesh::Mesh Mesh;

        //- Type of the internal field from which this GeometricField is derived
        typedef DimensionedField<Type, GeoMesh> Internal;

        //- Type of the patch field of which the Boundary is composed
        typedef PatchField<Type> Patch;

        //- Type of the boundary field
        typedef GeometricBoundaryField<Type, PatchField, GeoMesh> Boundary;


private:

    // Private Data

        //- Current time index.
        //  Used to trigger the storing of the old-time value
        mutable label timeIndex_;

        //- Pointer to old time field
        mutable GeometricField<Type, PatchField, GeoMesh>* field0Ptr_;

        //-  Pointer to previous iteration (used for under-relaxation)
        mutable GeometricField<Type, PatchField, GeoMesh>* fieldPrevIterPtr_;

        //- Boundary Type field containing boundary field values
        Boundary boundaryField_;


    // Private Member Functions

        //- Read the field from the dictionary
        void readFields(const dictionary&);

        //- Read the field - create the field dictionary on-the-fly
        void readFields();

        //- Store the old-time field
        void storeOldTime() const;

        //- Find oldest time and set to nullObject
        void nullOldestTimeFound();


public:

    //- Runtime type information
    TypeName("GeometricField");


    // Public Typedefs

        typedef typename Field<Type>::cmptType cmptType;

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
            const word& patchFieldType=PatchField<Type>::calculatedType()
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
            const wordList& actualPatchTypes = wordList()
        );

        //- Constructor given IOobject, mesh, dimensioned<Type>
        //  and patch field type.
        GeometricField
        (
            const IOobject&,
            const Mesh&,
            const dimensioned<Type>&,
            const word& patchFieldType=PatchField<Type>::calculatedType()
        );

        //- Constructor given IOobject, mesh, dimensioned<Type>
        //  and patch field types.
        GeometricField
        (
            const IOobject&,
            const Mesh&,
            const dimensioned<Type>&,
            const wordList& wantedPatchTypes,
            const wordList& actualPatchTypes = wordList()
        );

        //- Constructor from components
        GeometricField
        (
            const IOobject&,
            const Internal&,
            const PtrList<PatchField<Type>>&
        );

        //- Constructor from components
        GeometricField
        (
            const IOobject&,
            const Mesh&,
            const dimensionSet&,
            const Field<Type>&,
            const PtrList<PatchField<Type>>&
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

        //- Construct as copy of tmp<GeometricField> deleting argument
        GeometricField(const tmp<GeometricField<Type, PatchField, GeoMesh>>&);

        //- Construct as copy resetting IO parameters
        GeometricField
        (
            const IOobject&,
            const GeometricField<Type, PatchField, GeoMesh>&
        );

        //- Construct as copy of tmp<GeometricField> resetting IO parameters
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

        //- Construct as copy resetting name
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
            const wordList& actualPatchTypes = wordList()
        );

        //- Construct as copy resetting IO parameters and boundary types
        GeometricField
        (
            const IOobject&,
            const tmp<GeometricField<Type, PatchField, GeoMesh>>&,
            const wordList& patchFieldTypes,
            const wordList& actualPatchTypes = wordList()
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
            const PtrList<PatchField<Type>>&
        );

        //- Return a temporary field constructed from name, mesh, dimensionSet
        //  and patch field type.
        static tmp<GeometricField<Type, PatchField, GeoMesh>> New
        (
            const word& name,
            const Mesh&,
            const dimensionSet&,
            const word& patchFieldType=PatchField<Type>::calculatedType()
        );

        //- Return a temporary field constructed from mesh, dimensioned<Type>
        //  and patch field type.
        static tmp<GeometricField<Type, PatchField, GeoMesh>> New
        (
            const word& name,
            const Mesh&,
            const dimensioned<Type>&,
            const word& patchFieldType=PatchField<Type>::calculatedType()
        );

        //- Return a temporary field constructed from mesh, dimensioned<Type>
        //  and patch field types.
        static tmp<GeometricField<Type, PatchField, GeoMesh>> New
        (
            const word& name,
            const Mesh&,
            const dimensioned<Type>&,
            const wordList& patchFieldTypes,
            const wordList& actualPatchTypes = wordList()
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
            const word&
        );

        //- Rename and reset patch fields types of temporary field and return
        static tmp<GeometricField<Type, PatchField, GeoMesh>> New
        (
            const word& newName,
            const tmp<GeometricField<Type, PatchField, GeoMesh>>&,
            const wordList& patchFieldTypes,
            const wordList& actualPatchTypes = wordList()
        );


    //- Destructor
    virtual ~GeometricField();


    // Member Functions

        //- Return a reference to the dimensioned internal field
        //  Note: this increments the event counter and checks the
        //  old-time fields; avoid in loops.
        Internal& ref();

        //- Return a const-reference to the dimensioned internal field
        inline const Internal& internalField() const;

        //- Return a const-reference to the dimensioned internal field
        //  of a "vol" field.  Useful in the formulation of source-terms
        //  for FV equations
        inline const Internal& v() const;

        //- Return a reference to the internal field
        //  Note: this increments the event counter and checks the
        //  old-time fields; avoid in loops.
        typename Internal::FieldType& primitiveFieldRef();

        //- Return a const-reference to the  internal field
        inline const typename Internal::FieldType& primitiveField() const;

        //- Return a reference to the boundary field
        //  Note: this increments the event counter and checks the
        //  old-time fields; avoid in loops.
        Boundary& boundaryFieldRef();

        //- Return const-reference to the boundary field
        inline const Boundary& boundaryField() const;

        //- Return the time index of the field
        inline label timeIndex() const;

        //- Return the time index of the field
        inline label& timeIndex();

        //- Return whether or not this is an old-time field
        bool isOldTime() const;

        //- Store the old-time fields
        void storeOldTimes() const;

        //- Return the number of old time fields stored
        label nOldTimes() const;

        //- Return old time field
        const GeometricField<Type, PatchField, GeoMesh>& oldTime() const;

        //- Return non-const old time field
        //  (Not a good idea but it is used for sub-cycling)
        GeometricField<Type, PatchField, GeoMesh>& oldTime();

        //- Return the n-th old time field
        const GeometricField<Type, PatchField, GeoMesh>& oldTime
        (
            const label n
        ) const;

        //- Return the n-th non-const old time field
        //  (Not a good idea but it is used for sub-cycling)
        GeometricField<Type, PatchField, GeoMesh>& oldTime(const label n);

        //- Store the field as the previous iteration value
        void storePrevIter() const;

        //- Return previous iteration field
        const GeometricField<Type, PatchField, GeoMesh>& prevIter() const;

        //- Delete old time and previous iteration fields
        void clearOldTimes();

        //- Set oldest time field pointer to nullObjectPtr
        //  to remove the oldest time field while maintaining a tag
        //  so that it is reinstated on the next call to oldTime().
        //
        //  This is used to temporarily remove the oldest-time field
        //  prior to mesh topology change.
        void nullOldestTime();

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
