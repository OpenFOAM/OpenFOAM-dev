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

#include "LagrangianMesh.H"
#include "DimensionedField.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Field>
inline void Foam::LagrangianMesh::checkFieldSize(const Field& field) const
{
    if (field.size() != size())
    {
        FatalErrorInFunction
            << "size of " << field.objectPath() << " = "
            << field.size() << " is not the same as the size of the "
            << "mesh = " << size() << abort(FatalError);
    }
}


template<class Type>
inline void Foam::LagrangianMesh::checkPtr
(
    const autoPtr<Type>& ptr,
    const word& name
) const
{
    if (!ptr.valid())
    {
        FatalErrorInFunction
            << name << " requested for a Lagrangian mesh in which "
            << "they have not been constructed"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline const Foam::polyMesh& Foam::LagrangianMesh::mesh() const
{
    return mesh_;
}


inline const Foam::barycentricIODynamicField&
Foam::LagrangianMesh::coordinates() const
{
    return coordinates_;
}


inline const Foam::labelIODynamicField& Foam::LagrangianMesh::celli() const
{
    return celli_;
}


inline const Foam::labelIODynamicField& Foam::LagrangianMesh::facei() const
{
    return facei_;
}


inline const Foam::labelIODynamicField& Foam::LagrangianMesh::faceTrii() const
{
    return faceTrii_;
}


inline const Foam::LagrangianBoundaryMesh&
Foam::LagrangianMesh::boundary() const
{
    return boundary_;
}


inline bool Foam::LagrangianMesh::changing() const
{
    return statesPtr_.valid();
}


inline const Foam::List<Foam::LagrangianState>&
Foam::LagrangianMesh::states() const
{
    checkPtr(statesPtr_, "States");
    return statesPtr_();
}


inline Foam::List<Foam::LagrangianState>& Foam::LagrangianMesh::states()
{
    checkPtr(statesPtr_, "States");
    return statesPtr_();
}


inline Foam::LagrangianState Foam::LagrangianMesh::state(const label i) const
{
    return statesPtr_.valid() ? statesPtr_()[i] : LagrangianState::none;
}


inline Foam::LagrangianState Foam::LagrangianMesh::state
(
    const LagrangianSubMesh& subMesh,
    const label subi
) const
{
    return state(subMesh.start() + subi);
}


inline Foam::label Foam::LagrangianMesh::subMeshIndex() const
{
    if (statesPtr_.valid())
    {
        subMeshIndex_ ++;
    }
    else
    {
        subMeshIndex_ --;
    }

    if (subMeshIndex_ == labelMax || subMeshIndex_ == labelMin)
    {
        FatalErrorInFunction
            << "Sub mesh index counter has overflowed."
            << exit(FatalError);
    }

    return subMeshIndex_;
}


inline Foam::label Foam::LagrangianMesh::size() const
{
    return coordinates().size();
}


inline Foam::label Foam::LagrangianMesh::globalSize() const
{
    return returnReduce(size(), sumOp<label>());
}


inline Foam::label Foam::LagrangianMesh::size(const LagrangianMesh& mesh)
{
    return mesh.coordinates().size();
}


inline Foam::label Foam::LagrangianMesh::nStates() const
{
    return
        static_cast<label>(LagrangianState::onPatchZero)
      + mesh().boundaryMesh().size() + 1;
}


inline Foam::label Foam::LagrangianMesh::nGroups() const
{
    return
        static_cast<label>(LagrangianGroup::onPatchZero)
      + mesh().boundaryMesh().size() + 1;
}


inline Foam::label Foam::LagrangianMesh::stateToGroupi
(
    const LagrangianState state
) const
{
    return
        min
        (
            static_cast<label>(stateToGroup(state)),
            nGroups() - 1
        );
}


inline Foam::LagrangianSubMesh Foam::LagrangianMesh::sub
(
    const LagrangianGroup group
) const
{
    checkPtr(offsetsPtr_, "Offsets");
    return LagrangianSubMesh(*this, offsetsPtr_(), group);
}


inline Foam::LagrangianSubMesh Foam::LagrangianMesh::subNone() const
{
    return LagrangianSubMesh(*this, LagrangianGroup::none, 0, size());
}


inline const Foam::LagrangianSubMesh& Foam::LagrangianMesh::subAll() const
{
    return subAll_;
}


inline Foam::LagrangianSubMesh Foam::LagrangianMesh::subIncomplete() const
{
    checkPtr(offsetsPtr_, "Offsets");
    return
        LagrangianSubMesh
        (
            *this,
            LagrangianGroup::none,
            size()
          - offsetsPtr_()[static_cast<label>(LagrangianGroup::inInternalMesh)],
            offsetsPtr_()[static_cast<label>(LagrangianGroup::inInternalMesh)]
        );
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline const Foam::vector&
Foam::LagrangianMesh::linearDisplacement::operator()(const label i) const
{
    return linear_[i];
}


inline const Foam::vector& Foam::LagrangianMesh::linearDisplacement::operator()
(
    const label i,
    const scalar f
) const
{
    return linear_[i];
}


inline Foam::Pair<Foam::vector>
Foam::LagrangianMesh::parabolicDisplacement::operator()(const label i) const
{
    return Pair<vector>(linear_[i], quadratic_[i]);
}


inline Foam::vector Foam::LagrangianMesh::parabolicDisplacement::operator()
(
    const label i,
    const scalar f
) const
{
    return linear_[i] + (1 - f)*quadratic_[i];
}


// ************************************************************************* //
