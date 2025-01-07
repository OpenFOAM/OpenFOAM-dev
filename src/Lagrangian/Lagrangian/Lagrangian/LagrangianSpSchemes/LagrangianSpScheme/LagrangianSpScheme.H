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
    Foam::Lagrangian::SpScheme

Description
    Abstract base class for Lagrangian Sp schemes

SourceFiles
    LagrangianSpScheme.C
    LagrangianSpSchemes.C

\*---------------------------------------------------------------------------*/

#ifndef LagrangianSpScheme_H
#define LagrangianSpScheme_H

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
                          Class SpScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class SpType>
class SpScheme
:
    public tmp<SpScheme<Type, SpType>>::refCount
{
protected:

    // Protected Data

        //- Reference to the mesh
        const LagrangianMesh& mesh_;


    // Protected Member Functions

        //- Generate the inner product of a scalar coefficient with a field
        static tmp<LagrangianSubField<Type>> inner
        (
            const LagrangianSubScalarField& Sp,
            const LagrangianSubSubField<Type>& psi
        );

        //- Generate the inner product of a tensor coefficient with a field
        static tmp
        <
            LagrangianSubField
            <
                typename innerProduct<tensor, Type>::type
            >
        > inner
        (
            const LagrangianSubTensorField& Sp,
            const LagrangianSubSubField<Type>& psi
        );


public:

    //- Runtime type information
    TypeName("SpScheme");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        tmp,
        SpScheme,
        Istream,
        (const LagrangianMesh& mesh, Istream& is),
        (mesh, is)
    );


    // Constructors

        //- Construct from a mesh
        SpScheme(const LagrangianMesh& mesh)
        :
            mesh_(mesh)
        {}

        //- Construct from a mesh and a stream
        SpScheme(const LagrangianMesh& mesh, Istream&)
        :
            mesh_(mesh)
        {}

        //- Disallow default bitwise copy construction
        SpScheme(const SpScheme&) = delete;


    // Selectors

        //- Return a pointer to a new SpScheme
        static tmp<SpScheme<Type, SpType>> New
        (
            const LagrangianMesh& mesh,
            Istream& is
        );


    //- Destructor
    virtual ~SpScheme();


    // Member Functions

        //- Return mesh reference
        const LagrangianMesh& mesh() const
        {
            return mesh_;
        }

        //- Return the source matrix
        virtual tmp<LagrangianEqn<Type>> LagrangianmSp
        (
            const LagrangianSubField<SpType>& Sp,
            const LagrangianSubSubField<Type>& psi
        ) = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const SpScheme&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Lagrangian
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineLagrangianSpScheme(Type, SpType)                                 \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        namespace Lagrangian                                                   \
        {                                                                      \
            typedef SpScheme<Type, SpType> Type##SpType##SpScheme;             \
                                                                               \
            defineNamedTemplateTypeNameAndDebug                                \
            (                                                                  \
                Type##SpType##SpScheme,                                        \
                0                                                              \
            );                                                                 \
                                                                               \
            defineTemplateRunTimeSelectionTable                                \
            (                                                                  \
                Type##SpType##SpScheme,                                        \
                Istream                                                        \
            );                                                                 \
        }                                                                      \
    }


#define makeLagrangianSpScheme(Type, SpType, SpSchemeType)                     \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        namespace Lagrangian                                                   \
        {                                                                      \
            typedef SpScheme<Type, SpType> Type##SpType##SpScheme;             \
                                                                               \
            namespace SpSchemes                                                \
            {                                                                  \
                typedef SpSchemeType<Type, SpType> Type##SpType##SpSchemeType; \
                                                                               \
                defineNamedTemplateTypeNameAndDebug                            \
                (                                                              \
                    Type##SpType##SpSchemeType,                                \
                    0                                                          \
                );                                                             \
                                                                               \
                addToRunTimeSelectionTable                                     \
                (                                                              \
                    Type##SpType##SpScheme,                                    \
                    Type##SpType##SpSchemeType,                                \
                    Istream                                                    \
                );                                                             \
            }                                                                  \
        }                                                                      \
    }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "LagrangianSpScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
