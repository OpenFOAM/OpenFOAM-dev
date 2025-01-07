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
    Foam::Lagrangian::accumulationScheme

Description
    Abstract base class for Lagrangian accumulation schemes

SourceFiles
    LagrangianAccumulationScheme.C
    LagrangianAccumulationSchemes.C

\*---------------------------------------------------------------------------*/

#ifndef LagrangianAccumulationScheme_H
#define LagrangianAccumulationScheme_H

#include "tmp.H"
#include "dimensionedType.H"
#include "LagrangianFieldsFwd.H"
#include "LagrangianSubFieldsFwd.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class LagrangianMesh;

template<class Type>
class LagrangianEqn;

namespace Lagrangian
{

/*---------------------------------------------------------------------------*\
                          Class accumulationScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class accumulationScheme
:
    public tmp<accumulationScheme<Type>>::refCount
{
private:

    // Private Member Functions

        //- Accumulate the Lagrangian field into a primitive cell field
        virtual void accumulate
        (
            const LagrangianSubSubField<Type>& lPsi,
            Field<Type>& cPsi
        ) = 0;


protected:

    // Protected Data

        //- Reference to the mesh
        const LagrangianMesh& mesh_;


public:

    //- Runtime type information
    TypeName("accumulationScheme");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        tmp,
        accumulationScheme,
        Istream,
        (const LagrangianMesh& mesh, Istream& is),
        (mesh, is)
    );


    // Constructors

        //- Construct from a mesh
        accumulationScheme(const LagrangianMesh& mesh)
        :
            mesh_(mesh)
        {}

        //- Construct from a mesh and a stream
        accumulationScheme(const LagrangianMesh& mesh, Istream&)
        :
            mesh_(mesh)
        {}

        //- Disallow default bitwise copy construction
        accumulationScheme(const accumulationScheme&) = delete;


    // Selectors

        //- Return a pointer to a new accumulationScheme
        static tmp<accumulationScheme<Type>> New
        (
            const LagrangianMesh& mesh,
            Istream& is
        );


    //- Destructor
    virtual ~accumulationScheme();


    // Member Functions

        //- Return mesh reference
        const LagrangianMesh& mesh() const
        {
            return mesh_;
        }

        //- Accumulate a Lagrangian field into a finite volume field
        template<class CellMesh, template<class> class PrimitiveField>
        tmp<DimensionedField<Type, CellMesh>> accumulate
        (
            const DimensionedField<Type, LagrangianMesh, PrimitiveField>& lPsi
        );

        //- Accumulate a Lagrangian sub-field into a finite volume field
        template<class CellMesh, template<class> class PrimitiveField>
        void accumulate
        (
            const LagrangianSubField<Type, PrimitiveField>& lPsi,
            DimensionedField<Type, CellMesh>& vPsi
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const accumulationScheme&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Lagrangian
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineLagrangianAccumulationScheme(Type, nullArg)                      \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        namespace Lagrangian                                                   \
        {                                                                      \
            typedef accumulationScheme<Type> Type##AccumulationScheme;         \
                                                                               \
            defineNamedTemplateTypeNameAndDebug                                \
            (                                                                  \
                Type##AccumulationScheme,                                      \
                0                                                              \
            );                                                                 \
                                                                               \
            defineTemplateRunTimeSelectionTable                                \
            (                                                                  \
                Type##AccumulationScheme,                                      \
                Istream                                                        \
            );                                                                 \
        }                                                                      \
    }


#define makeLagrangianAccumulationScheme(Type, AccumulationSchemeType)         \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        namespace Lagrangian                                                   \
        {                                                                      \
            typedef accumulationScheme<Type> Type##AccumulationScheme;         \
                                                                               \
            namespace accumulationSchemes                                      \
            {                                                                  \
                typedef AccumulationSchemeType<Type>                           \
                    Type##AccumulationSchemeType;                              \
                                                                               \
                defineNamedTemplateTypeNameAndDebug                            \
                (                                                              \
                    Type##AccumulationSchemeType,                              \
                    0                                                          \
                );                                                             \
                                                                               \
                addToRunTimeSelectionTable                                     \
                (                                                              \
                    Type##AccumulationScheme,                                  \
                    Type##AccumulationSchemeType,                              \
                    Istream                                                    \
                );                                                             \
            }                                                                  \
        }                                                                      \
    }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "LagrangianAccumulationScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
