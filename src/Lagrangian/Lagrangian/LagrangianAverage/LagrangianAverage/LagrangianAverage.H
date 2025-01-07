/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::LagrangianAverage

Description
    Base class for methods of local-averaging of Lagrangian properties

SourceFiles
    LagrangianAverage.C
    LagrangianAverages.C

\*---------------------------------------------------------------------------*/

#ifndef LagrangianAverage_H
#define LagrangianAverage_H

#include "dimensionSet.H"
#include "LagrangianFieldsFwd.H"
#include "LagrangianSubFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
class LagrangianAverage;

template<class CellMesh, class Type>
struct LagrangianAverageNewReturnType
{
    typedef autoPtr<LagrangianAverage<Type>> type;
};

template<class Type>
struct LagrangianAverageNewReturnType<LagrangianMesh, Type>
{
    // No definition. SFINAE then resolves the New templates below.
};

/*---------------------------------------------------------------------------*\
                      Class LagrangianAverage Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class LagrangianAverage
{
    // Private Data

        //- The name of the result
        const word name_;

        //- Reference to the mesh
        const LagrangianMesh& mesh_;

        //- Dimensions of the result
        const dimensionSet dimensions_;


    // Private Member Functions

        //- Interpolate into a sub-field, where possible. Calling code deals
        //  with default values where no elements are available to interpolate.
        virtual void interpolate(LagrangianSubField<Type>& result) const = 0;


public:

    //- Runtime type information
    TypeName("LagrangianAverage");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        LagrangianAverage,
        dictionary,
        (
            const word& name,
            const LagrangianMesh& mesh,
            const dimensionSet& dimensions,
            const Field<scalar>& cellWeightSum
        ),
        (name, mesh, dimensions, cellWeightSum)
    );


    // Constructors

        //- Construct with a name, for a mesh and with given dimensions
        LagrangianAverage
        (
            const word& name,
            const LagrangianMesh& mesh,
            const dimensionSet& dimensions
        );


    // Selectors

        //- Select to average a field
        template<template<class> class WeightPF, template<class> class PsiPF>
        static autoPtr<LagrangianAverage<Type>> New
        (
            const word& type,
            const word& name,
            const DimensionedField<scalar, LagrangianMesh, WeightPF>& weight,
            const DimensionedField<Type, LagrangianMesh, PsiPF>& psi
        );

        //- Select to average a field with a given cell weight sum
        template<class CellMesh, template<class> class WeightPsiPF>
        static typename LagrangianAverageNewReturnType<CellMesh, Type>::type New
        (
            const word& type,
            const word& name,
            const DimensionedField<scalar, CellMesh>& cellWeightSum,
            const DimensionedField<Type, LagrangianMesh, WeightPsiPF>& weightPsi
        );


    //- Destructor
    virtual ~LagrangianAverage();


    // Member Functions

        //- Remove weighted values from the average
        virtual void remove
        (
            const LagrangianSubSubField<scalar>& weight,
            const LagrangianSubSubField<Type>& psi
        ) = 0;

        //- Remove values from the average
        void remove(const LagrangianSubSubField<Type>& weightPsi);

        //- Add weighted values to the average
        virtual void add
        (
            const LagrangianSubSubField<scalar>& weight,
            const LagrangianSubSubField<Type>& psi,
            const bool cache
        ) = 0;

        //- Add values to the average
        void add
        (
            const LagrangianSubSubField<Type>& weightPsi,
            const bool cache
        );

        //- Correct weighted values in the average
        virtual void correct
        (
            const LagrangianSubSubField<scalar>& weight,
            const LagrangianSubSubField<Type>& psi,
            const bool cache
        ) = 0;

        //- Correct values in the average
        void correct
        (
            const LagrangianSubSubField<Type>& weightPsi,
            const bool cache
        );

        //- Interpolate to a sub-mesh
        tmp<LagrangianSubField<Type>> interpolate
        (
            const LagrangianSubMesh& subMesh
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const LagrangianAverage<Type>&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "LagrangianAverage.C"
#endif


#define defineLagrangianAverage(Type, nullArg)                                 \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        typedef LagrangianAverage<Type> Type##Average;                         \
                                                                               \
        defineNamedTemplateTypeNameAndDebug                                    \
        (                                                                      \
            Type##Average,                                                     \
            0                                                                  \
        );                                                                     \
                                                                               \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            Type##Average,                                                     \
            dictionary                                                         \
        );                                                                     \
    }


#define makeLagrangianAverage(Type, AverageType)                               \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        typedef LagrangianAverage<Type> Type##Average;                         \
                                                                               \
        namespace LagrangianAverages                                           \
        {                                                                      \
            typedef AverageType<Type> Type##AverageType;                       \
                                                                               \
            defineNamedTemplateTypeNameAndDebug                                \
            (                                                                  \
                Type##AverageType,                                             \
                0                                                              \
            );                                                                 \
                                                                               \
            addToRunTimeSelectionTable                                         \
            (                                                                  \
                Type##Average,                                                 \
                Type##AverageType,                                             \
                dictionary                                                     \
            );                                                                 \
        }                                                                      \
    }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
