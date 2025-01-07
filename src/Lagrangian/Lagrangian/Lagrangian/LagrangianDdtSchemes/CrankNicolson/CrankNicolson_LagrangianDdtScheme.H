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
    Foam::Lagrangian::ddtSchemes::CrankNicolson

Description
    Basic first-order implicit CrankNicolson Lagrangian ddt scheme

SourceFiles
    CrankNicolson_LagrangianDdtScheme.C
    CrankNicolson_LagrangianDdtSchemes.C

\*---------------------------------------------------------------------------*/

#ifndef CrankNicolson_LagrangianDdtScheme_H
#define CrankNicolson_LagrangianDdtScheme_H

#include "LagrangianDdtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Lagrangian
{
namespace ddtSchemes
{

/*---------------------------------------------------------------------------*\
                          Class CrankNicolson Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class CrankNicolson
:
    public ddtScheme<Type>
{
    // Private Member Functions

        //- Return the name of the deltaTSp0 field for the given field
        word deltaTSp0Name(const LagrangianSubSubField<Type>& psi) const;

        //- Return the name of the S0 field for the given field
        word S0Name(const LagrangianSubSubField<Type>& psi) const;


public:

    //- Runtime type information
    TypeName("CrankNicolson");


    // Constructors

        //- Inherit constructors
        using ddtScheme<Type>::ddtScheme;


    // Member Functions

        //- Initialise time-derivative information
        virtual bool LagrangianmInitDdt
        (
            const dimensionSet& mDims,
            const LagrangianSubSubField<Type>& psi,
            const bool instantaneousDdt
        );

        //- Return the no-time-derivative matrix
        virtual tmp<LagrangianEqn<Type>> LagrangianmNoDdt
        (
            const LagrangianSubScalarField& deltaT,
            const dimensionSet& mDims,
            const LagrangianSubSubField<Type>& psi
        );

        //- Return the time-derivative matrix
        virtual tmp<LagrangianEqn<Type>> LagrangianmDdt
        (
            const LagrangianSubScalarField& deltaT,
            LagrangianSubSubField<Type>& psi
        );

        //- Return the time-derivative matrix
        virtual tmp<LagrangianEqn<Type>> LagrangianmDdt
        (
            const LagrangianSubScalarField& deltaT,
            const LagrangianSubScalarSubField& m,
            LagrangianSubSubField<Type>& psi
        );

        //- Return the instantaneous time-derivative
        virtual tmp<LagrangianSubField<Type>> LagrangiancDdt
        (
            const LagrangianSubSubField<Type>& psi
        );

        //- Return the instantaneous time-derivative
        virtual tmp<LagrangianSubField<Type>> LagrangiancDdt
        (
            const LagrangianSubScalarSubField& m,
            const LagrangianSubSubField<Type>& psi
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace ddtSchemes
} // End namespace Lagrangian
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CrankNicolson_LagrangianDdtScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
