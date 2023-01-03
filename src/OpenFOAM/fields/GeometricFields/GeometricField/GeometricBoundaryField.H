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
    Foam::GeometricBoundaryField

Description
    Generic GeometricBoundaryField class.

SourceFiles
    GeometricBoundaryField.C

\*---------------------------------------------------------------------------*/

#ifndef GeometricBoundaryField_H
#define GeometricBoundaryField_H

#include "dimensionedTypes.H"
#include "DimensionedField.H"
#include "FieldField.H"
#include "lduInterfaceFieldPtrsList.H"
#include "LduInterfaceFieldPtrsList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class dictionary;

/*---------------------------------------------------------------------------*\
                       Class GeometricBoundaryField Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class GeoMesh, template<class> class PrimitiveField>
class GeometricBoundaryField
:
    public FieldField<GeoMesh::template PatchField, Type>
{
public:

    // Public Typedefs

        //- Type of boundary mesh on which this boundary is instantiated
        typedef typename GeoMesh::BoundaryMesh BoundaryMesh;

        //- Type of the internal field from which this GeometricField is derived
        typedef DimensionedField<Type, GeoMesh, PrimitiveField> Internal;

        //- Type of the patch field of which the Boundary is composed
        typedef typename GeoMesh::template PatchField<Type> Patch;


private:

    // Private Data

        //- Reference to BoundaryMesh for which this field is defined
        const BoundaryMesh& bmesh_;


public:

    //- Declare friendship with other geometric boundary fields
    template<class Type2, class GeoMesh2, template<class> class PrimitiveField2>
    friend class GeometricBoundaryField;


    // Constructors

        //- Construct from a BoundaryMesh
        GeometricBoundaryField(const BoundaryMesh&);

        //- Construct from a BoundaryMesh,
        //  reference to the internal field
        //  and a patch field type
        GeometricBoundaryField
        (
            const BoundaryMesh&,
            const Internal&,
            const word&
        );

        //- Construct from a BoundaryMesh,
        //  reference to the internal field
        //  and a wordList of patch field types and optionally
        //  the actual patch types (to override constraint patches)
        GeometricBoundaryField
        (
            const BoundaryMesh&,
            const Internal&,
            const wordList& wantedPatchTypes,
            const wordList& actualPatchTypes = wordList()
        );

        //- Construct from a BoundaryMesh,
        //  reference to the internal field
        //  and a PtrList<Patch>
        GeometricBoundaryField
        (
            const BoundaryMesh&,
            const Internal&,
            const PtrList<Patch>&
        );

        //- Construct as copy setting the reference to the internal field
        GeometricBoundaryField(const Internal&, const GeometricBoundaryField&);

        //- Construct as copy setting the reference to the internal field
        template<template<class> class PrimitiveField2>
        explicit GeometricBoundaryField
        (
            const Internal& field,
            const GeometricBoundaryField<Type, GeoMesh, PrimitiveField2>& btf
        );

        //- Copy constructor deleted
        //  as it would not set the internalField reference correctly
        GeometricBoundaryField(const GeometricBoundaryField&) = delete;

        //- Move constructor deleted
        //  as it would not set the internalField reference correctly
        GeometricBoundaryField(GeometricBoundaryField&&) = delete;

        //- Construct from dictionary
        GeometricBoundaryField
        (
            const BoundaryMesh&,
            const Internal&,
            const dictionary&
        );


    // Member Functions

        //- Read the boundary field
        void readField(const Internal& field, const dictionary& dict);

        //- Update the boundary condition coefficients
        void updateCoeffs();

        //- Evaluate boundary conditions
        void evaluate();

        //- Return a list of the patch field types
        wordList types() const;

        //- Return BoundaryField of the cell values neighbouring
        //  the boundary
        tmp<GeometricBoundaryField> boundaryInternalField() const;

        //- Return BoundaryField of the values on the other side of couples
        tmp<GeometricBoundaryField> boundaryNeighbourField() const;

        //- Return a list of pointers for each patch field with only those
        //  pointing to interfaces being set
        LduInterfaceFieldPtrsList<Type> interfaces() const;

        //- Return a list of pointers for each patch field with only those
        //  pointing to interfaces being set
        lduInterfaceFieldPtrsList scalarInterfaces() const;

        //- Reset the boundary field contents to the given field
        //  Used for mesh to mesh mapping
        void reset(const GeometricBoundaryField&);

        //- Write boundary field as dictionary entry
        void writeEntry(const word& keyword, Ostream& os) const;


    // Member Operators

        //- Assignment operator
        void operator=(const GeometricBoundaryField&);

        //- Move assignment operator
        void operator=(GeometricBoundaryField&&);

        //- Assignment to FieldField<GeoMesh::template PatchField, Type>
        void operator=(const FieldField<GeoMesh::template PatchField, Type>&);

        //- Assignment to FieldField<OtherPatchField, Type>
        template<template<class> class OtherPatchField>
        void operator=(const FieldField<OtherPatchField, Type>&);

        //- Assignment to Type
        void operator=(const Type&);

        //- Forced assignment to
        //  BoundaryField<Type, GeoMesh::template PatchField, BoundaryMesh>
        void operator==(const GeometricBoundaryField&);

        //- Forced assignment to FieldField<GeoMesh::template PatchField, Type>
        void operator==(const FieldField<GeoMesh::template PatchField, Type>&);

        //- Forced assignment to FieldField<OtherPatchField, Type>
        template<template<class> class OtherPatchField>
        void operator==(const FieldField<OtherPatchField, Type>&);

        //- Forced assignment to Type
        void operator==(const Type&);
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

#ifdef NoRepository
    #include "GeometricBoundaryField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
