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
    Foam::DimensionedFunction2

Description
    Wrapper around a Function2 which provides evaluation for dimensioned and
    geometric fields

SourceFiles
    DimensionedFunction2.C

\*---------------------------------------------------------------------------*/

#ifndef DimensionedFunction2_H
#define DimensionedFunction2_H

#include "GeometricField.H"
#include "Function2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class DimensionedFunction2 Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class DimensionedFunction2
{
    // Private Data

        //- The dimensions of the x-coordinate
        const dimensionSet xDimensions_;

        //- The dimensions of the y-coordinate
        const dimensionSet yDimensions_;

        //- The dimensions of the value
        const dimensionSet valueDimensions_;

        //- The function
        const autoPtr<Function2<Type>> function_;


    // Private Constructors

        //- Construct from a name, units and a dictionary
        DimensionedFunction2
        (
            const word& name,
            const Function2s::unitSets& units,
            const dictionary& dict
        );


    // Private Member Functions

        //- Set a dimensioned field to the value
        template
        <
            class GeoMesh,
            template<class> class ResultPrimitiveField,
            template<class> class XPrimitiveField,
            template<class> class YPrimitiveField
        >
        void value
        (
            DimensionedField<scalar, GeoMesh, ResultPrimitiveField>& result,
            const DimensionedField<scalar, GeoMesh, XPrimitiveField>& x,
            const DimensionedField<scalar, GeoMesh, YPrimitiveField>& y
        ) const;

        //- Set a geometric field to the value
        template
        <
            class GeoMesh,
            template<class> class ResultPrimitiveField,
            template<class> class XPrimitiveField,
            template<class> class YPrimitiveField
        >
        void value
        (
            GeometricField<scalar, GeoMesh, ResultPrimitiveField>& result,
            const GeometricField<scalar, GeoMesh, XPrimitiveField>& x,
            const GeometricField<scalar, GeoMesh, YPrimitiveField>& y
        ) const;


public:

    // Selectors

        //- Select from dictionary
        static autoPtr<DimensionedFunction2<Type>> New
        (
            const word& name,
            const Function2s::unitSets& units,
            const dictionary& dict
        );

        //- Select from dictionary
        static autoPtr<DimensionedFunction2<Type>> New
        (
            const word& name,
            const unitSet& xUnits,
            const unitSet& yUnits,
            const unitSet& valueUnits,
            const dictionary& dict
        );


    // Member Functions

        //- Return the value as a field
        tmp<Field<Type>> value
        (
            const scalarField& x,
            const scalarField& y
        ) const;

        //- Return the value as a dimensioned field
        template
        <
            class GeoMesh,
            template<class> class XPrimitiveField,
            template<class> class YPrimitiveField
        >
        tmp<DimensionedField<Type, GeoMesh>> value
        (
            const DimensionedField<scalar, GeoMesh, XPrimitiveField>& x,
            const DimensionedField<scalar, GeoMesh, YPrimitiveField>& y
        ) const;

        //- Return the value as a geometric field
        template
        <
            class GeoMesh,
            template<class> class XPrimitiveField,
            template<class> class YPrimitiveField
        >
        tmp<GeometricField<Type, GeoMesh>> value
        (
            const GeometricField<scalar, GeoMesh, XPrimitiveField>& x,
            const GeometricField<scalar, GeoMesh, YPrimitiveField>& y
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DimensionedFunction2.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
