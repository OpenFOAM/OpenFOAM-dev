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
    Foam::Lagrangian::SpSchemes::implicit

Description
    Basic implicit Lagrangian Sp scheme

SourceFiles
    implicit_LagrangianSpScheme.C
    implicit_LagrangianSpSchemes.C

\*---------------------------------------------------------------------------*/

#ifndef implicit_LagrangianSpScheme_H
#define implicit_LagrangianSpScheme_H

#include "LagrangianSpScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Lagrangian
{
namespace SpSchemes
{

/*---------------------------------------------------------------------------*\
                          Class implicit Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class SpType>
class implicit
:
    public SpScheme<Type, SpType>
{
public:

    //- Runtime type information
    TypeName("implicit");


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
    #include "implicit_LagrangianSpScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
