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

InNamespace
    Foam::Lagrangianm

Description
    Functions for calculating sources in a Lagrangian equation

SourceFiles
    LagrangianmSp.C

\*---------------------------------------------------------------------------*/

#ifndef LagrangianmSp_H
#define LagrangianmSp_H

#include "LagrangianFieldsFwd.H"
#include "LagrangianSubFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
class LagrangianEqn;

/*---------------------------------------------------------------------------*\
                    Namespace Lagrangianm functions Declaration
\*---------------------------------------------------------------------------*/

namespace Lagrangianm
{
    template<class Type, class SpType>
    tmp<LagrangianEqn<Type>> Sp
    (
        const LagrangianSubField<SpType>& Sp,
        const LagrangianSubSubField<Type>& psi
    );

    template<class Type, class SpType, template<class> class PrimitiveField>
    tmp<LagrangianEqn<Type>> Sp
    (
        const LagrangianSubField<SpType>& Sp,
        const LagrangianSubField<Type, PrimitiveField>& psi
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "LagrangianmSp.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
