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

#include "LagrangiancAccumulate.H"
#include "LagrangianAccumulationScheme.H"
#include "LagrangianMesh.H"
#include "LagrangianSubFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class CellMesh, class Type, template<class> class PrimitiveField>
Foam::tmp<Foam::DimensionedField<Type, CellMesh>>
Foam::Lagrangianc::accumulate
(
    const DimensionedField<Type, LagrangianMesh, PrimitiveField>& lPsi
)
{
    return accumulate<CellMesh>(lPsi, lPsi.name());
}


template<class CellMesh, class Type, template<class> class PrimitiveField>
Foam::tmp<Foam::DimensionedField<Type, CellMesh>>
Foam::Lagrangianc::accumulate
(
    const DimensionedField<Type, LagrangianMesh, PrimitiveField>& lPsi,
    const word& name
)
{
    const LagrangianMesh& mesh = lPsi.mesh();

    return
        Lagrangian::accumulationScheme<Type>::New
        (
            mesh,
            mesh.schemes().accumulation(name)
        ).ref().template accumulate<CellMesh>(lPsi);
}


template<class CellMesh, class Type, template<class> class PrimitiveField>
void Foam::Lagrangianc::accumulate
(
    const LagrangianSubField<Type, PrimitiveField>& lPsi,
    DimensionedField<Type, CellMesh>& vPsi
)
{
    accumulate<CellMesh>(lPsi, vPsi, lPsi.name());
}


template<class CellMesh, class Type, template<class> class PrimitiveField>
void Foam::Lagrangianc::accumulate
(
    const LagrangianSubField<Type, PrimitiveField>& lPsi,
    DimensionedField<Type, CellMesh>& vPsi,
    const word& name
)
{
    const LagrangianMesh& mesh = lPsi.mesh().mesh();

    Lagrangian::accumulationScheme<Type>::New
    (
        mesh,
        mesh.schemes().accumulation(name)
    ).ref().template accumulate<CellMesh>(lPsi, vPsi);
}


// ************************************************************************* //
