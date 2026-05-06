/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::FunctionalGeometricField

Description
    GeometricField with a corresponding run-time selected function
    to evaluate and update the field.

See also
    Foam::GeometricField
    Foam::DimensionedFieldFunction
    Foam::UFunctionalDimensionedField

SourceFiles
    FunctionalGeometricField.C

\*---------------------------------------------------------------------------*/

#ifndef FunctionalGeometricField_H
#define FunctionalGeometricField_H

#include "FunctionalGeometricFieldFwd.H"
#include "GeometricField.H"
#include "DimensionedFieldFunction.H"
#include "SlicedDimensionedField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions
template<class Type, class GeoMesh>
void writeEntry
(
    Ostream& os,
    const FunctionalGeometricField<Type, GeoMesh>& udff
);

/*---------------------------------------------------------------------------*\
                    Class FunctionalGeometricField Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class GeoMesh>
class FunctionalGeometricField
:
    public GeometricField<Type, GeoMesh>
{
    // Private Type Definitions

        //- The type of the patch
        typedef typename GeoMesh::BoundaryMesh::value_type GeoPatch;


    // Private Member Functions

        //- Set the patch field pointers
        void setPatchFields();

        //- Read the functions
        bool readFuncs(const dictionary& dict);

        //- Clone the functions
        void cloneFuncs(const FunctionalGeometricField<Type, GeoMesh>& udff);


    // Private Data

        //- Name of the function entry/sub dictionary
        const word funcName_;

        //- Optional internal function pointer
        autoPtr<DimensionedFieldFunction<DimensionedField<Type, GeoMesh>>>
            internalFuncPtr_;

        //- Optional patch field pointers
        PtrList<SlicedDimensionedField<Type, GeoPatch>>
            patchFieldPtrs_;

        //- Optional patch function pointers
        PtrList<DimensionedFieldFunction<DimensionedField<Type, GeoPatch>>>
            patchFuncPtrs_;

        //- Optional default value
        const Type defaultValue_;


public:

    // Constructors

        //- Construct from name, mesh, dimensions, field and dictionary
        FunctionalGeometricField
        (
            const word& name,
            const word& funcName,
            const GeoMesh&,
            const dimensionSet& dimensions,
            const dictionary&
        );

        //- Construct from name, mesh, dimensions, field, dictionary
        //  and default value
        FunctionalGeometricField
        (
            const word& name,
            const word& funcName,
            const GeoMesh&,
            const dimensionSet& dimensions,
            const dictionary&,
            const Type& defaultValue
        );

        //- Construct from given FunctionalGeometricField
        //  with new mesh and field
        FunctionalGeometricField
        (
            const FunctionalGeometricField<Type, GeoMesh>&,
            const GeoMesh&
        );

        //- Copy constructor
        FunctionalGeometricField
        (
            const FunctionalGeometricField<Type, GeoMesh>&
        );


    // Member Functions

        //- Reset
        void reset();

        //- Update
        bool update();

        using regIOobject::write;

        //- Write
        void write(Ostream&) const;


    // Global functions

        friend void writeEntry <Type, GeoMesh>
        (
            Ostream& os,
            const FunctionalGeometricField<Type, GeoMesh>& udff
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "FunctionalGeometricField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
