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
#include "Explicit_LagrangianSpScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class SpType>
Foam::tmp<Foam::LagrangianEqn<Type>> Foam::Lagrangianm::Sp
(
    const LagrangianSubField<SpType>& Sp,
    const LagrangianSubSubField<Type>& psi
)
{
    if (isNull(Sp))
    {
        return tmp<LagrangianEqn<Type>>(new LagrangianEqn<Type>(psi));
    }

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


template<class Type>
Foam::tmp<Foam::LagrangianEqn<Type>> Foam::Lagrangianm::explicitSp0
(
    const LagrangianSubScalarField& Sp,
    const LagrangianSubSubField<Type>& psi
)
{
    if (isNull(Sp))
    {
        return tmp<LagrangianEqn<Type>>(new LagrangianEqn<Type>(psi));
    }

    tmp<LagrangianSubSubField<Type>> tpsi0 =
        psi.nOldTimes()
      ? toSubField(psi.oldTime())
      : tmp<LagrangianSubSubField<Type>>(psi);

    return
        Lagrangian::SpSchemes::Explicit<Type, scalar>
        (
            psi.mesh().mesh()
        ).LagrangianmSp(Sp, tpsi0());
}


template<class Type>
Foam::tmp<Foam::LagrangianEqn<Type>> Foam::Lagrangianm::explicitSp0
(
    const LagrangianSubScalarField& Sp,
    const LagrangianSubScalarSubField& m,
    const LagrangianSubSubField<Type>& psi
)
{
    if (isNull(Sp))
    {
        return tmp<LagrangianEqn<Type>>(new LagrangianEqn<Type>(psi));
    }

    tmp<LagrangianSubScalarSubField> tm0 =
        m.nOldTimes()
      ? toSubField(m.oldTime())
      : tmp<LagrangianSubSubField<Type>>(m);
    tmp<LagrangianSubSubField<Type>> tpsi0 =
        psi.nOldTimes()
      ? toSubField(psi.oldTime())
      : tmp<LagrangianSubSubField<Type>>(psi);

    return
        Lagrangian::SpSchemes::Explicit<Type, scalar>
        (
            psi.mesh().mesh()
        ).LagrangianmSp(Sp*tm0, tpsi0());
}


template<class Type>
Foam::tmp<Foam::LagrangianEqn<Type>> Foam::Lagrangianm::implicitDeltaTSp
(
    const LagrangianSubScalarField& deltaTSp,
    const LagrangianSubSubField<Type>& psi
)
{
    if (isNull(deltaTSp))
    {
        return tmp<LagrangianEqn<Type>>(new LagrangianEqn<Type>(psi));
    }

    tmp<LagrangianEqn<Type>> tEqn(new LagrangianEqn<Type>(psi));

    tEqn.ref().deltaTSp += deltaTSp;

    return tEqn;
}


// ************************************************************************* //
