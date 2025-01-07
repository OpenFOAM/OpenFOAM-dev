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
    Foam::Lagrangian::ddtScheme

Description
    Abstract base class for Lagrangian ddt schemes

SourceFiles
    LagrangianDdtScheme.C
    LagrangianDdtSchemes.C

\*---------------------------------------------------------------------------*/

#ifndef LagrangianDdtScheme_H
#define LagrangianDdtScheme_H

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
                          Class ddtScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class ddtScheme
:
    public tmp<ddtScheme<Type>>::refCount
{
protected:

    // Protected Data

        //- Reference to the mesh
        const LagrangianMesh& mesh_;


public:

    //- Runtime type information
    TypeName("ddtScheme");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        tmp,
        ddtScheme,
        Istream,
        (const LagrangianMesh& mesh, Istream& is),
        (mesh, is)
    );


    // Constructors

        //- Construct from a mesh
        ddtScheme(const LagrangianMesh& mesh)
        :
            mesh_(mesh)
        {}

        //- Construct from a mesh and a stream
        ddtScheme(const LagrangianMesh& mesh, Istream&)
        :
            mesh_(mesh)
        {}

        //- Disallow default bitwise copy construction
        ddtScheme(const ddtScheme&) = delete;


    // Selectors

        //- Return a pointer to a new ddtScheme
        static tmp<ddtScheme<Type>> New
        (
            const LagrangianMesh& mesh,
            Istream& is
        );


    //- Destructor
    virtual ~ddtScheme();


    // Member Functions

        //- Return mesh reference
        const LagrangianMesh& mesh() const
        {
            return mesh_;
        }

        //- Initialise time-derivative information
        virtual bool LagrangianmInitDdt
        (
            const dimensionSet& mDims,
            const LagrangianSubSubField<Type>& psi,
            const bool instantaneousDdt
        ) = 0;

        //- Return the no-time-derivative matrix
        virtual tmp<LagrangianEqn<Type>> LagrangianmNoDdt
        (
            const LagrangianSubScalarField& deltaT,
            const dimensionSet& mDims,
            const LagrangianSubSubField<Type>& psi
        ) = 0;

        //- Return the time-derivative matrix
        virtual tmp<LagrangianEqn<Type>> LagrangianmDdt
        (
            const LagrangianSubScalarField& deltaT,
            LagrangianSubSubField<Type>& psi
        ) = 0;

        //- Return the time-derivative matrix
        virtual tmp<LagrangianEqn<Type>> LagrangianmDdt
        (
            const LagrangianSubScalarField& deltaT,
            const LagrangianSubScalarSubField& m,
            LagrangianSubSubField<Type>& psi
        ) = 0;

        //- Return the instantaneous time-derivative
        virtual tmp<LagrangianSubField<Type>> LagrangiancDdt
        (
            const LagrangianSubSubField<Type>& psi
        ) = 0;

        //- Return the instantaneous time-derivative
        virtual tmp<LagrangianSubField<Type>> LagrangiancDdt
        (
            const LagrangianSubScalarSubField& m,
            const LagrangianSubSubField<Type>& psi
        ) = 0;

        //- Return the time-derivative matrix
        static tmp<LagrangianEqn<Type>> Lagrangianmddt
        (
            const LagrangianSubScalarField& deltaT,
            const LagrangianSubSubField<Type>& psi
        );

        //- Return the time-derivative matrix
        static tmp<LagrangianEqn<Type>> Lagrangianmddt
        (
            const LagrangianSubScalarField& deltaT,
            const LagrangianSubScalarSubField& m,
            const LagrangianSubSubField<Type>& psi
        );

        //- Return the explicit/forward time-derivative matrix
        static tmp<LagrangianEqn<Type>> Lagrangianmddt0
        (
            const LagrangianSubScalarField& deltaT,
            const LagrangianSubSubField<Type>& psi
        );

        //- Return the explicit/forward time-derivative matrix
        static tmp<LagrangianEqn<Type>> Lagrangianmddt0
        (
            const LagrangianSubScalarField& deltaT,
            const LagrangianSubScalarSubField& m,
            const LagrangianSubSubField<Type>& psi
        );


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ddtScheme&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Lagrangian
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineLagrangianDdtScheme(Type, nullArg)                               \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        namespace Lagrangian                                                   \
        {                                                                      \
            typedef ddtScheme<Type> Type##DdtScheme;                           \
                                                                               \
            defineNamedTemplateTypeNameAndDebug(Type##DdtScheme, 0);           \
                                                                               \
            defineTemplateRunTimeSelectionTable(Type##DdtScheme, Istream);     \
        }                                                                      \
    }


#define makeLagrangianDdtScheme(Type, DdtSchemeType)                           \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        namespace Lagrangian                                                   \
        {                                                                      \
            typedef ddtScheme<Type> Type##DdtScheme;                           \
                                                                               \
            namespace ddtSchemes                                               \
            {                                                                  \
                typedef DdtSchemeType<Type> Type##DdtSchemeType;               \
                                                                               \
                defineNamedTemplateTypeNameAndDebug(Type##DdtSchemeType, 0);   \
                                                                               \
                addToRunTimeSelectionTable                                     \
                (                                                              \
                    Type##DdtScheme,                                           \
                    Type##DdtSchemeType,                                       \
                    Istream                                                    \
                );                                                             \
            }                                                                  \
        }                                                                      \
    }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "LagrangianDdtScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
