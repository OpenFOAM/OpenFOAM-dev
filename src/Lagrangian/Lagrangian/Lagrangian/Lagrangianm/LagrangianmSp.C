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

#include "LagrangianmSp.H"
#include "LagrangianSpScheme.H"
#include "LagrangianMesh.H"
#include "LagrangianSubFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class SpType>
Foam::tmp<Foam::LagrangianEqn<Type>> Foam::Lagrangianm::Sp
(
    const LagrangianSubField<SpType>& Sp,
    const LagrangianSubSubField<Type>& psi
)
{
    const LagrangianMesh& mesh = psi.mesh().mesh();

    return
        Lagrangian::SpScheme<Type, SpType>::New
        (
            mesh,
            mesh.schemes().Sp("Sp(" + Sp.name() + ',' + psi.name() + ')')
        ).ref().LagrangianmSp(Sp, psi);
}


template<class Type, class SpType, template<class> class PrimitiveField>
Foam::tmp<Foam::LagrangianEqn<Type>> Foam::Lagrangianm::Sp
(
    const LagrangianSubField<SpType>& Sp,
    const LagrangianSubField<Type, PrimitiveField>& psi
)
{
    return Lagrangianm::Sp(Sp, toSubField(psi)());
}


// ************************************************************************* //
