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
    Foam::Lagrangian::SpSchemes::Explicit

Description
    Basic explicit Lagrangian Sp scheme

SourceFiles
    Explicit_LagrangianSpScheme.C
    Explicit_LagrangianSpSchemes.C

\*---------------------------------------------------------------------------*/

#ifndef Explicit_LagrangianSpScheme_H
#define Explicit_LagrangianSpScheme_H

#include "LagrangianSpScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Lagrangian
{
namespace SpSchemes
{

/*---------------------------------------------------------------------------*\
                          Class Explicit Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class SpType>
class Explicit
:
    public SpScheme<Type, SpType>
{
public:

    //- Runtime type information
    TypeName("explicit");


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
    #include "Explicit_LagrangianSpScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
