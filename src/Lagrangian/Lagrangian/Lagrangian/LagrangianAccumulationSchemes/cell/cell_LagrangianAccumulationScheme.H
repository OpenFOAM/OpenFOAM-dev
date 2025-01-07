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
    Basic cell Lagrangian accumulation scheme

SourceFiles
    cell_LagrangianAccumulationScheme.C
    cell_LagrangianAccumulationSchemes.C

\*---------------------------------------------------------------------------*/

#ifndef cell_LagrangianAccumulationScheme_H
#define cell_LagrangianAccumulationScheme_H

#include "LagrangianAccumulationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Lagrangian
{
namespace accumulationSchemes
{

/*---------------------------------------------------------------------------*\
                          Class cell Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class cell
:
    public accumulationScheme<Type>
{
private:

    // Private Member Functions

        //- Accumulate the Lagrangian field into a primitive cell field
        virtual void accumulate
        (
            const LagrangianSubSubField<Type>& lPsi,
            Field<Type>& cPsi
        );


public:

    //- Runtime type information
    TypeName("cell");


    // Constructors

        //- Inherit constructors
        using accumulationScheme<Type>::accumulationScheme;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace accumulationSchemes
} // End namespace Lagrangian
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "cell_LagrangianAccumulationScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
