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

#include "CloudStateField.H"
#include "toSubField.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
Foam::LagrangianSubSubField<Type>& Foam::CloudStateField<Type>::ref
(
    const LagrangianSubMesh& subMesh
)
{
    // Evaluate and store if it doesn't already exist for the sub-mesh
    if (!psiSubSubPtr_.valid() || psiSubSubMeshIndex_ != subMesh.index())
    {
        psiSubSubPtr_.reset(subMesh.sub(*this).ptr());

        psiSubSubMeshIndex_ = subMesh.index();
    }

    return psiSubSubPtr_();
}


template<class Type>
void Foam::CloudStateField<Type>::clear()
{
    psiSubSubPtr_.clear();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
const Foam::LagrangianDynamicField<Type>&
Foam::CloudStateField<Type>::operator()(const LagrangianMesh&) const
{
    return *this;
}


template<class Type>
Foam::tmp<Foam::LagrangianSubSubField<Type>>
Foam::CloudStateField<Type>::operator()
(
    const LagrangianSubMesh& subMesh
) const
{
    return subMesh.sub(*this);
}


template<class Type>
Foam::tmp<Foam::LagrangianSubField<Type>>
Foam::CloudStateField<Type>::operator()
(
    const LagrangianModel& model,
    const LagrangianSubMesh& subMesh
) const
{
    return this->sources()[model.name()].value(model, subMesh);
}


template<class Type>
Foam::tmp<Foam::LagrangianSubSubField<Type>>
Foam::CloudStateField<Type>::operator()
(
    const LagrangianModelRef& model,
    const LagrangianSubMesh& subMesh
) const
{
    if (!model.valid())
    {
        return subMesh.sub(*this);
    }
    else
    {
        return toSubField(operator()(model(), subMesh));
    }
}


// ************************************************************************* //
