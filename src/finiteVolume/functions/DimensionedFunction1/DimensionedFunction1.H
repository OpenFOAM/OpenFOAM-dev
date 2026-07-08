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
    Foam::DimensionedFunction1

Description
    Wrapper around a Function1 which provides evaluation for dimensioned and
    geometric fields

SourceFiles
    DimensionedFunction1.C

\*---------------------------------------------------------------------------*/

#ifndef DimensionedFunction1_H
#define DimensionedFunction1_H

#include "GeometricField.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class DimensionedFunction1 Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class DimensionedFunction1
{
    // Private Data

        //- The dimensions of the x-coordinate
        const dimensionSet xDimensions_;

        //- The dimensions of the value
        const dimensionSet valueDimensions_;

        //- The function
        const autoPtr<Function1<Type>> function_;


    // Private Constructors

        //- Construct from a name, units and a dictionary
        DimensionedFunction1
        (
            const word& name,
            const Function1s::unitSets& units,
            const dictionary& dict
        );


    // Private Member Functions

        //- Set a dimensioned field to the value
        template
        <
            class GeoMesh,
            template<class> class ResultPrimitiveField,
            template<class> class XPrimitiveField
        >
        void value
        (
            DimensionedField<scalar, GeoMesh, ResultPrimitiveField>& result,
            const DimensionedField<scalar, GeoMesh, XPrimitiveField>& x
        ) const;

        //- Set a dimensioned field to the derivative
        template
        <
            class GeoMesh,
            template<class> class ResultPrimitiveField,
            template<class> class XPrimitiveField
        >
        void derivative
        (
            DimensionedField<scalar, GeoMesh, ResultPrimitiveField>& result,
            const DimensionedField<scalar, GeoMesh, XPrimitiveField>& x
        ) const;

        //- Set a geometric field to the value
        template
        <
            class GeoMesh,
            template<class> class ResultPrimitiveField,
            template<class> class XPrimitiveField
        >
        void value
        (
            GeometricField<scalar, GeoMesh, ResultPrimitiveField>& result,
            const GeometricField<scalar, GeoMesh, XPrimitiveField>& x
        ) const;

        //- Set a geometric field to the derivative
        template
        <
            class GeoMesh,
            template<class> class ResultPrimitiveField,
            template<class> class XPrimitiveField
        >
        void derivative
        (
            GeometricField<scalar, GeoMesh, ResultPrimitiveField>& result,
            const GeometricField<scalar, GeoMesh, XPrimitiveField>& x
        ) const;


public:

    // Selectors

        //- Select from dictionary
        static autoPtr<DimensionedFunction1<Type>> New
        (
            const word& name,
            const Function1s::unitSets& units,
            const dictionary& dict
        );

        //- Select from dictionary
        static autoPtr<DimensionedFunction1<Type>> New
        (
            const word& name,
            const unitSet& xUnits,
            const unitSet& valueUnits,
            const dictionary& dict
        );


    // Member Functions

        //- Return the value as a field
        tmp<Field<Type>> value(const scalarField& x) const;

        //- Return the derivative as a field
        tmp<Field<Type>> derivative(const scalarField& x) const;

        //- Return the value as a dimensioned field
        template<class GeoMesh, template<class> class XPrimitiveField>
        tmp<DimensionedField<Type, GeoMesh>> value
        (
            const DimensionedField<scalar, GeoMesh, XPrimitiveField>& x
        ) const;

        //- Return the derivative as a dimensioned field
        template<class GeoMesh, template<class> class XPrimitiveField>
        tmp<DimensionedField<Type, GeoMesh>> derivative
        (
            const DimensionedField<scalar, GeoMesh, XPrimitiveField>& x
        ) const;

        //- Return the value as a dimensioned field
        template<class GeoMesh, template<class> class XPrimitiveField>
        tmp<DimensionedField<Type, GeoMesh>> value
        (
            const tmp<DimensionedField<scalar, GeoMesh, XPrimitiveField>>& tx
        ) const;

        //- Return the derivative as a dimensioned field
        template<class GeoMesh, template<class> class XPrimitiveField>
        tmp<DimensionedField<Type, GeoMesh>> derivative
        (
            const tmp<DimensionedField<scalar, GeoMesh, XPrimitiveField>>& tx
        ) const;

        //- Return the value as a geometric field
        template<class GeoMesh, template<class> class XPrimitiveField>
        tmp<GeometricField<Type, GeoMesh>> value
        (
            const GeometricField<scalar, GeoMesh, XPrimitiveField>& x
        ) const;

        //- Return the derivative as a geometric field
        template<class GeoMesh, template<class> class XPrimitiveField>
        tmp<GeometricField<Type, GeoMesh>> derivative
        (
            const GeometricField<scalar, GeoMesh, XPrimitiveField>& x
        ) const;

        //- Return the value as a geometric field
        template<class GeoMesh, template<class> class XPrimitiveField>
        tmp<GeometricField<Type, GeoMesh>> value
        (
            const tmp<GeometricField<scalar, GeoMesh, XPrimitiveField>>& tx
        ) const;

        //- Return the derivative as a geometric field
        template<class GeoMesh, template<class> class XPrimitiveField>
        tmp<GeometricField<Type, GeoMesh>> derivative
        (
            const tmp<GeometricField<scalar, GeoMesh, XPrimitiveField>>& tx
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DimensionedFunction1.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
