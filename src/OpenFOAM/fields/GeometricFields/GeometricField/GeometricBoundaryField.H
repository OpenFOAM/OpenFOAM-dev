/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class dictionary;

// Forward declaration of friend functions and operators

template<class Type, template<class> class PatchField, class GeoMesh>
class GeometricField;

/*---------------------------------------------------------------------------*\
                       Class GeometricBoundaryField Declaration
\*---------------------------------------------------------------------------*/

template<class Type, template<class> class PatchField, class GeoMesh>
class GeometricBoundaryField
:
    public FieldField<PatchField, Type>
{
public:

    // Public Typedefs

        //- Type of boundary mesh on which this boundary is instantiated
        typedef typename GeoMesh::BoundaryMesh BoundaryMesh;

        //- Type of the internal field from which this GeometricField is derived
        typedef DimensionedField<Type, GeoMesh> Internal;


private:

    // Private Data

        //- Reference to BoundaryMesh for which this field is defined
        const BoundaryMesh& bmesh_;


public:

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
        //  and a PtrList<PatchField<Type>>
        GeometricBoundaryField
        (
            const BoundaryMesh&,
            const Internal&,
            const PtrList<PatchField<Type>>&
        );

        //- Construct as copy setting the reference to the internal field
        GeometricBoundaryField(const Internal&, const GeometricBoundaryField&);

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
        void reset(const GeometricBoundaryField<Type, PatchField, GeoMesh>&);

        //- Write boundary field as dictionary entry
        void writeEntry(const word& keyword, Ostream& os) const;


    // Member Operators

        //- Assignment operator
        void operator=(const GeometricBoundaryField&);

        //- Move assignment operator
        void operator=(GeometricBoundaryField&&);

        //- Assignment to FieldField<PatchField, Type>
        void operator=(const FieldField<PatchField, Type>&);

        //- Assignment to FieldField<OtherPatchField, Type>
        template<template<class> class OtherPatchField>
        void operator=(const FieldField<OtherPatchField, Type>&);

        //- Assignment to Type
        void operator=(const Type&);


        //- Forced assignment to
        //  BoundaryField<Type, PatchField, BoundaryMesh>
        void operator==(const GeometricBoundaryField&);

        //- Forced assignment to FieldField<PatchField, Type>
        void operator==(const FieldField<PatchField, Type>&);

        //- Forced assignment to FieldField<OtherPatchField, Type>
        template<template<class> class OtherPatchField>
        void operator==(const FieldField<OtherPatchField, Type>&);

        //- Forced assignment to Type
        void operator==(const Type&);
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

#ifdef NoRepository
    #include "GeometricBoundaryField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
