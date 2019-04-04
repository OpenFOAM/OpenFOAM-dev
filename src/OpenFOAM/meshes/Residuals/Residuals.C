/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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

#include "Residuals.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Residuals<Type>::Residuals(const polyMesh& mesh)
:
    MeshObject<polyMesh, GeometricMeshObject, Residuals<Type>>(mesh),
    prevTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
Foam::List<Foam::word> Foam::Residuals<Type>::fieldNames(const polyMesh& mesh)
{
    return MeshObject<polyMesh, GeometricMeshObject, Residuals<Type>>::New
    (
        mesh
    ).HashTable<DynamicList<SolverPerformance<Type>>>::toc();
}


template<class Type>
bool Foam::Residuals<Type>::found(const polyMesh& mesh, const word& fieldName)
{
    return MeshObject<polyMesh, GeometricMeshObject, Residuals<Type>>::New
    (
        mesh
    ).HashTable<DynamicList<SolverPerformance<Type>>>::found(fieldName);
}


template<class Type>
const Foam::DynamicList<Foam::SolverPerformance<Type>>&
Foam::Residuals<Type>::field
(
    const polyMesh& mesh,
    const word& fieldName
)
{
    return MeshObject<polyMesh, GeometricMeshObject, Residuals<Type>>::New
    (
        mesh
    )[fieldName];
}


template<class Type>
void Foam::Residuals<Type>::append
(
    const polyMesh& mesh,
    const SolverPerformance<Type>& sp
)
{
    Residuals<Type>& residuals = const_cast<Residuals<Type>&>
    (
        MeshObject<polyMesh, GeometricMeshObject, Residuals<Type>>::New
        (
            mesh
        )
    );

    HashTable<DynamicList<SolverPerformance<Type>>>& table = residuals;

    const label timeIndex =
        mesh.time().subCycling()
      ? mesh.time().prevTimeState().timeIndex()
      : mesh.time().timeIndex();

    if (residuals.prevTimeIndex_ != timeIndex)
    {
        // Reset solver performance between iterations
        residuals.prevTimeIndex_ = timeIndex;
        table.clear();
    }

    if (table.found(sp.fieldName()))
    {
        table[sp.fieldName()].append(sp);
    }
    else
    {
        table.insert
        (
            sp.fieldName(),
            DynamicList<SolverPerformance<Type>>(1, sp)
        );
    }
}


// ************************************************************************* //
