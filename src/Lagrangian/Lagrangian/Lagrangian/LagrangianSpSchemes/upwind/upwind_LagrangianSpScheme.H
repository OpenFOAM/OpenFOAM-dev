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
    Foam::Lagrangian::SpSchemes::upwind

Description
    Upwind/Su-Sp Lagrangian Sp scheme

SourceFiles
    upwind_LagrangianSpScheme.C
    upwind_LagrangianSpSchemes.C

\*---------------------------------------------------------------------------*/

#ifndef upwind_LagrangianSpScheme_H
#define upwind_LagrangianSpScheme_H

#include "LagrangianSpScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Lagrangian
{
namespace SpSchemes
{

/*---------------------------------------------------------------------------*\
                          Class upwind Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class SpType>
class upwind
:
    public SpScheme<Type, SpType>
{
public:

    //- Runtime type information
    TypeName("upwind");


    // Constructors

        //- Inherit constructors
        using SpScheme<Type, SpType>::SpScheme;


    // Member Functions

        //- Return the source matrix
        virtual tmp<LagrangianEqn<Type>> LagrangianmSp
        (
            const LagrangianSubField<SpType>& Sp,
            const LagrangianSubSubField<Type>& psi
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace SpSchemes
} // End namespace Lagrangian
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "upwind_LagrangianSpScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
