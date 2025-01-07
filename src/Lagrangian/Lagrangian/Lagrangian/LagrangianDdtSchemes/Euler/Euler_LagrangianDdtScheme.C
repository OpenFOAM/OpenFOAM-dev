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

\*---------------------------------------------------------------------------*/

#include "Euler_LagrangianDdtScheme.H"
#include "LagrangianFields.H"
#include "LagrangianEqn.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
bool Foam::Lagrangian::ddtSchemes::Euler<Type>::LagrangianmInitDdt
(
    const dimensionSet& mDims,
    const LagrangianSubSubField<Type>& psi,
    const bool instantaneousDdt
)
{
    if (instantaneousDdt)
    {
        FatalErrorInFunction
            << typeName << " " << ddtScheme<Type>::typeName << " does not "
            << "support the generation of the instantaneous time derivative"
            << exit(FatalError);
    }

    return false;
}


template<class Type>
Foam::tmp<Foam::LagrangianEqn<Type>>
Foam::Lagrangian::ddtSchemes::Euler<Type>::LagrangianmNoDdt
(
    const LagrangianSubScalarField& deltaT,
    const dimensionSet& mDims,
    const LagrangianSubSubField<Type>& psi
)
{
    return tmp<LagrangianEqn<Type>>(new LagrangianEqn<Type>(deltaT, psi));
}


template<class Type>
Foam::tmp<Foam::LagrangianEqn<Type>>
Foam::Lagrangian::ddtSchemes::Euler<Type>::LagrangianmDdt
(
    const LagrangianSubScalarField& deltaT,
    LagrangianSubSubField<Type>& psi
)
{
    tmp<LagrangianEqn<Type>> tEqn =
        Lagrangian::ddtScheme<Type>::Lagrangianmddt(deltaT, psi);

    // Set the equation name and the non-const field reference
    tEqn->op(LagrangianEqn<Type>(psi));

    return tEqn;
}


template<class Type>
Foam::tmp<Foam::LagrangianEqn<Type>>
Foam::Lagrangian::ddtSchemes::Euler<Type>::LagrangianmDdt
(
    const LagrangianSubScalarField& deltaT,
    const LagrangianSubScalarSubField& m,
    LagrangianSubSubField<Type>& psi
)
{
    tmp<LagrangianEqn<Type>> tEqn =
        Lagrangian::ddtScheme<Type>::Lagrangianmddt(deltaT, m, psi);

    // Set the equation name and the non-const field reference
    tEqn->op(LagrangianEqn<Type>(psi));

    return tEqn;
}


template<class Type>
Foam::tmp<Foam::LagrangianSubField<Type>>
Foam::Lagrangian::ddtSchemes::Euler<Type>::LagrangiancDdt
(
    const LagrangianSubSubField<Type>& psi
)
{
    LagrangianmInitDdt(dimless, psi, true);
    return tmp<LagrangianSubField<Type>>(nullptr);
}


template<class Type>
Foam::tmp<Foam::LagrangianSubField<Type>>
Foam::Lagrangian::ddtSchemes::Euler<Type>::LagrangiancDdt
(
    const LagrangianSubScalarSubField& m,
    const LagrangianSubSubField<Type>& psi
)
{
    LagrangianmInitDdt(dimless, psi, true);
    return tmp<LagrangianSubField<Type>>(nullptr);
}


// ************************************************************************* //
