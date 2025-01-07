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
    Foam::Lagrangian::ddtSchemes::Euler

Description
    Basic first-order implicit Euler Lagrangian ddt scheme

SourceFiles
    Euler_LagrangianDdtScheme.C
    Euler_LagrangianDdtSchemes.C

\*---------------------------------------------------------------------------*/

#ifndef Euler_LagrangianDdtScheme_H
#define Euler_LagrangianDdtScheme_H

#include "LagrangianDdtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Lagrangian
{
namespace ddtSchemes
{

/*---------------------------------------------------------------------------*\
                          Class Euler Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Euler
:
    public ddtScheme<Type>
{
public:

    //- Runtime type information
    TypeName("Euler");


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
    #include "Euler_LagrangianDdtScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
