/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025-2026 OpenFOAM Foundation
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

#include "CloudStateFieldRef.H"
#include "toSubField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::CloudStateFieldRef<Type>::CloudStateFieldRef
(
    LagrangianDynamicField<Type>& ref
)
:
    ref_(ref),
    psiAllPtr_(ref_.mesh().subAll().sub(ref_).ptr())
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
const Foam::word& Foam::CloudStateFieldRef<Type>::name() const
{
    return ref_.name();
}


template<class Type>
const Foam::LagrangianMesh& Foam::CloudStateFieldRef<Type>::mesh() const
{
    return ref_.mesh();
}


template<class Type>
const Foam::dimensionSet& Foam::CloudStateFieldRef<Type>::dimensions() const
{
    return ref_.dimensions();
}


template<class Type>
const Foam::LagrangianSubSubField<Type>& Foam::CloudStateFieldRef<Type>::ref
(
    const LagrangianSubMesh& subMesh
) const
{
    // Error if this is the all-mesh
    if (&subMesh == &subMesh.mesh().subAll())
    {
        FatalErrorInFunction
            << "Non-constant access is not provided to the all-mesh sub-field"
            << exit(FatalError);
    }

    // Evaluate and store if it doesn't already exist for the sub-mesh
    if (!psiSubSubPtr_.valid() || psiSubSubMeshIndex_ != subMesh.index())
    {
        psiSubSubPtr_.reset(subMesh.sub(ref_).ptr());

        psiSubSubMeshIndex_ = subMesh.index();
    }

    return psiSubSubPtr_();
}


template<class Type>
Foam::LagrangianSubSubField<Type>& Foam::CloudStateFieldRef<Type>::ref
(
    const LagrangianSubMesh& subMesh
)
{
    // Evaluate and store if it doesn't already exist for the sub-mesh
    static_cast<const CloudStateFieldRef<Type>&>(*this).ref(subMesh);

    return psiSubSubPtr_();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
Foam::CloudStateFieldRef<Type>::operator
const Foam::LagrangianDynamicField<Type>&() const
{
    return ref_;
}


template<class Type>
Foam::CloudStateFieldRef<Type>::operator
LagrangianDynamicField<Type>&()
{
    return ref_;
}


template<class Type>
const Foam::LagrangianDynamicField<Type>&
Foam::CloudStateFieldRef<Type>::operator()(const LagrangianMesh&) const
{
    return ref_;
}


template<class Type>
const Foam::LagrangianSubSubField<Type>&
Foam::CloudStateFieldRef<Type>::operator()
(
    const LagrangianSubMesh& subMesh
) const
{
    // Update and return the all-mesh field
    if (&subMesh == &subMesh.mesh().subAll())
    {
        psiAllPtr_->UList<Type>::shallowCopy(ref_);

        return psiAllPtr_();
    }

    // Evaluate and store if it doesn't already exist for the sub-mesh
    if (!psiSubSubPtr_.valid() || psiSubSubMeshIndex_ != subMesh.index())
    {
        psiSubSubPtr_.reset(subMesh.sub(ref_).ptr());

        psiSubSubMeshIndex_ = subMesh.index();
    }

    return psiSubSubPtr_();
}


template<class Type>
Foam::tmp<Foam::LagrangianSubSubField<Type>>
Foam::CloudStateFieldRef<Type>::operator()
(
    const LagrangianModelRef& model,
    const LagrangianSubMesh& subMesh
) const
{
    if (model.isSource())
    {
        return
            toSubField
            (
                ref_.sources()[model.source().name()].value
                (
                    model.source(),
                    subMesh // !!! <-- Replace with model.S()
                )
            );
    }

    if (model.isInjection())
    {
        return
            toSubField
            (
                ref_.sources()[model.injection().name()].value
                (
                    model.injection(),
                    subMesh
                )
            );
    }

    return tmp<LagrangianSubSubField<Type>>(operator()(subMesh));
}


// ************************************************************************* //
