/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "regIOobject.H"
#include "Field.H"
#include "dimensionedType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Type, class GeoMesh> class DimensionedField;

template<class Type, class GeoMesh> Ostream& operator<<
(
    Ostream&,
    const DimensionedField<Type, GeoMesh>&
);

template<class Type, class GeoMesh> Ostream& operator<<
(
    Ostream&,
    const tmp<DimensionedField<Type, GeoMesh>>&
);


/*---------------------------------------------------------------------------*\
                      Class DimensionedField Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class GeoMesh>
class DimensionedField
:
    public regIOobject,
    public Field<Type>
{

public:

    // Public Typedefs

        //- Type of mesh on which this DimensionedField is instantiated
        typedef typename GeoMesh::Mesh Mesh;

        //- Type of the field from which this DimensionedField is derived
        typedef Field<Type> FieldType;

        //- Component type of the elements of the field
        typedef typename Field<Type>::cmptType cmptType;


private:

    // Private Data

        //- Reference to mesh
        const Mesh& mesh_;

        //- Dimension set for this field
        dimensionSet dimensions_;


    // Private Member Functions

        void readIfPresent(const word& fieldDictEntry = "value");


public:

    //- Runtime type information
    TypeName("DimensionedField");


    // Static Member Functions

        //- Return a null DimensionedField
        inline static const DimensionedField<Type, GeoMesh>& null();


    // Constructors

        //- Construct from components
        DimensionedField
        (
            const IOobject&,
            const Mesh& mesh,
            const dimensionSet&,
            const Field<Type>&
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
            const DimensionedField<Type, GeoMesh>&
        );

        //- Copy constructor or re-use as specified.
        DimensionedField
        (
            DimensionedField<Type, GeoMesh>&,
            bool reuse
        );

        //- Move constructor
        DimensionedField
        (
            DimensionedField<Type, GeoMesh>&&
        );

        //- Copy constructor of tmp<DimensionedField> deleting argument
        DimensionedField
        (
            const tmp<DimensionedField<Type, GeoMesh>>&
        );

        //- Copy constructor resetting IO parameters
        DimensionedField
        (
            const IOobject&,
            const DimensionedField<Type, GeoMesh>&
        );

        //- Copy constructor resetting IO parameters and re-use as specified.
        DimensionedField
        (
            const IOobject&,
            DimensionedField<Type, GeoMesh>&,
            bool reuse
        );

        //- Copy constructor resetting name
        DimensionedField
        (
            const word& newName,
            const DimensionedField<Type, GeoMesh>&
        );

        //- Copy constructor resetting name and re-use as specified.
        DimensionedField
        (
            const word& newName,
            DimensionedField<Type, GeoMesh>&,
            bool reuse
        );

        //- Copy constructor resetting name
        DimensionedField
        (
            const word& newName,
            const tmp<DimensionedField<Type, GeoMesh>>&
        );

        //- Clone
        tmp<DimensionedField<Type, GeoMesh>> clone() const;

        //- Return a temporary field constructed from name, mesh
        //  and dimensionSet
        static tmp<DimensionedField<Type, GeoMesh>> New
        (
            const word& name,
            const Mesh& mesh,
            const dimensionSet&
        );

        //- Return a temporary field constructed from name, mesh
        //  and dimensionType providing dimensions and value
        static tmp<DimensionedField<Type, GeoMesh>> New
        (
            const word& name,
            const Mesh& mesh,
            const dimensioned<Type>&
        );

        //- Return a temporary field constructed from name and a field
        static tmp<DimensionedField<Type, GeoMesh>> New
        (
            const word& newName,
            const DimensionedField<Type, GeoMesh>&
        );

        //- Return renamed temporary field
        static tmp<DimensionedField<Type, GeoMesh>> New
        (
            const word& newName,
            const tmp<DimensionedField<Type, GeoMesh>>&
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

        inline const Field<Type>& field() const;

        inline Field<Type>& field();

        //- Return a component field of the field
        tmp<DimensionedField<cmptType, GeoMesh>> component
        (
            const direction
        ) const;

        //- Replace a component field of the field
        void replace
        (
            const direction,
            const DimensionedField<cmptType, GeoMesh>&
        );

        //- Replace a component field of the field
        void replace
        (
            const direction,
            const tmp<DimensionedField<cmptType, GeoMesh>>&
        );

        //- Return the field transpose (only defined for second rank tensors)
        tmp<DimensionedField<Type, GeoMesh>> T() const;

        //- Calculate and return arithmetic average
        dimensioned<Type> average() const;

        //- Calculate and return weighted average
        dimensioned<Type> weightedAverage
        (
            const DimensionedField<scalar, GeoMesh>&
        ) const;

        //- Calculate and return weighted average
        dimensioned<Type> weightedAverage
        (
            const tmp<DimensionedField<scalar, GeoMesh>>&
        ) const;


        // Write

            bool writeData(Ostream&, const word& fieldDictEntry) const;

            bool writeData(Ostream&) const;


    // Member Operators

        void operator=(const DimensionedField<Type, GeoMesh>&);
        void operator=(DimensionedField<Type, GeoMesh>&&);
        void operator=(const tmp<DimensionedField<Type, GeoMesh>>&);
        void operator=(const dimensioned<Type>&);
        void operator=(const zero&);

        void operator+=(const DimensionedField<Type, GeoMesh>&);
        void operator+=(const tmp<DimensionedField<Type, GeoMesh>>&);

        void operator-=(const DimensionedField<Type, GeoMesh>&);
        void operator-=(const tmp<DimensionedField<Type, GeoMesh>>&);

        void operator*=(const DimensionedField<scalar, GeoMesh>&);
        void operator*=(const tmp<DimensionedField<scalar, GeoMesh>>&);

        void operator/=(const DimensionedField<scalar, GeoMesh>&);
        void operator/=(const tmp<DimensionedField<scalar, GeoMesh>>&);

        void operator+=(const dimensioned<Type>&);
        void operator-=(const dimensioned<Type>&);

        void operator*=(const dimensioned<scalar>&);
        void operator/=(const dimensioned<scalar>&);


    // Ostream Operators

        friend Ostream& operator<< <Type, GeoMesh>
        (
            Ostream&,
            const DimensionedField<Type, GeoMesh>&
        );

        friend Ostream& operator<< <Type, GeoMesh>
        (
            Ostream&,
            const tmp<DimensionedField<Type, GeoMesh>>&
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
