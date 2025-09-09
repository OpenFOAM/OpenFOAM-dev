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
    Functions for calculating the time derivative for a Lagrangian equation

SourceFiles
    LagrangianmDdt.C

\*---------------------------------------------------------------------------*/

#ifndef LagrangianmDdt_H
#define LagrangianmDdt_H

#include "LagrangianFieldsFwd.H"
#include "LagrangianSubFieldsFwd.H"
#include "dimensionSet.H"

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
    //- Initialise time-derivative information. Registers fields used to
    //  construct for higher-order schemes and to evaluate instantaneous
    //  time derivatives.
    template<class Type>
    bool initDdt
    (
        const dimensionSet& mDims,
        const LagrangianSubSubField<Type>& psi,
        const bool instantaneousDdt = false
    );

    //- Initialise time-derivative information. Registers fields used to
    //  construct for higher-order schemes and to evaluate instantaneous
    //  time derivatives.
    template<class Type, template<class> class PrimitiveField>
    bool initDdt
    (
        const dimensionSet& mDims,
        const LagrangianSubField<Type, PrimitiveField>& psi,
        const bool instantaneousDdt = false
    );

    //- Return the no-time-derivative matrix containing higher-order
    //  modifications to source terms
    template<class Type>
    tmp<LagrangianEqn<Type>> noDdt
    (
        const LagrangianSubScalarField& deltaT,
        const dimensionSet& mDims,
        const LagrangianSubSubField<Type>& psi
    );

    //- Return the no-time-derivative matrix containing higher-order
    //  modifications to source terms
    template<class Type, template<class> class PrimitiveField>
    tmp<LagrangianEqn<Type>> noDdt
    (
        const LagrangianSubScalarField& deltaT,
        const dimensionSet& mDims,
        const LagrangianSubField<Type, PrimitiveField>& psi
    );

    //- Return the time-derivative matrix
    template<class Type>
    tmp<LagrangianEqn<Type>> Ddt
    (
        const LagrangianSubScalarField& deltaT,
        LagrangianSubSubField<Type>& psi
    );

    //- Return the time-derivative matrix
    template<class Type>
    tmp<LagrangianEqn<Type>> Ddt
    (
        const LagrangianSubScalarField& deltaT,
        const LagrangianSubScalarSubField& m,
        LagrangianSubSubField<Type>& psi
    );

    //- Return a Euler implicit time derivative matrix
    template<class Type>
    tmp<LagrangianEqn<Type>> ddt
    (
        const LagrangianSubSubField<Type>& psi
    );

    //- Return a Euler implicit time derivative matrix
    template<class Type>
    tmp<LagrangianEqn<Type>> ddt
    (
        const LagrangianSubScalarSubField& m,
        const LagrangianSubSubField<Type>& psi
    );

    //- Return a Euler explicit/forward time derivative matrix
    template<class Type>
    tmp<LagrangianEqn<Type>> ddt0
    (
        const LagrangianSubSubField<Type>& psi
    );

    //- Return a Euler explicit/forward time derivative matrix
    template<class Type>
    tmp<LagrangianEqn<Type>> ddt0
    (
        const LagrangianSubScalarSubField& m,
        const LagrangianSubSubField<Type>& psi
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "LagrangianmDdt.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
