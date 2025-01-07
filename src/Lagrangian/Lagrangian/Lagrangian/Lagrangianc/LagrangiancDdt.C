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

#include "LagrangiancDdt.H"
#include "LagrangianDdtScheme.H"
#include "LagrangianMesh.H"
#include "LagrangianSubFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::LagrangianSubField<Type>> Foam::Lagrangianc::Ddt
(
    const LagrangianSubSubField<Type>& psi
)
{
    const LagrangianMesh& mesh = psi.mesh().mesh();

    return
        Lagrangian::ddtScheme<Type>::New
        (
            mesh,
            mesh.schemes().ddt("ddt(" + psi.name() + ')')
        ).ref().LagrangiancDdt(psi);
}


template<class Type>
Foam::tmp<Foam::LagrangianSubField<Type>> Foam::Lagrangianc::Ddt
(
    const LagrangianSubScalarSubField& m,
    const LagrangianSubSubField<Type>& psi
)
{
    const LagrangianMesh& mesh = psi.mesh().mesh();

    return
        Lagrangian::ddtScheme<Type>::New
        (
            mesh,
            mesh.schemes().ddt("ddt(" + m.name() + ',' + psi.name() + ')')
        ).ref().LagrangiancDdt(m, psi);
}


// ************************************************************************* //
