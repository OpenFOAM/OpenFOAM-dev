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

#ifndef LagrangianState_H
#define LagrangianState_H

#include "word.H"
#include "label.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Lagrangian state enumeration
enum class LagrangianState : label
{
    none = labelMin,
    complete = 0,
    inCell = 1,
    onInternalFace = 2,
    onPatchZero = 3,
    toBeRemoved = labelMax - 1
};

//- Return a string representation of a Lagrangian state enumeration
word name(const LagrangianState state);

//- Output operator for Lagrangian state
inline Ostream& operator<<(Ostream& os, const LagrangianState state)
{
    return os << name(state);
}

//- Lagrangian group enumeration
enum class LagrangianGroup : label
{
    none = labelMin,
    complete = 0,
    inInternalMesh = 1,
    onPatchZero = 2,
    toBeRemoved = labelMax - 1
};

//- Return a string representation of a Lagrangian group enumeration
word name(const LagrangianGroup group);

//- Output operator for Lagrangian group
inline Ostream& operator<<(Ostream& os, const LagrangianGroup group)
{
    return os << name(group);
}

//- Convert from a state enumeration to the corresponding group enumerations
inline LagrangianGroup stateToGroup(const LagrangianState& state)
{
    switch (state)
    {
        case LagrangianState::none:
            return LagrangianGroup::none;
        case LagrangianState::complete:
            return LagrangianGroup::complete;
        case LagrangianState::inCell:
            return LagrangianGroup::inInternalMesh;
        case LagrangianState::onInternalFace:
            return LagrangianGroup::inInternalMesh;
        default:
            return
                static_cast<LagrangianGroup>
                (
                    static_cast<label>(state)
                  - static_cast<label>(LagrangianState::onPatchZero)
                  + static_cast<label>(LagrangianGroup::onPatchZero)
                );
        case LagrangianState::toBeRemoved:
            return LagrangianGroup::toBeRemoved;
    }
}

//- Convert from a state enumeration to the corresponding group enumerations
inline LagrangianState groupToState(const LagrangianGroup& group)
{
    switch (group)
    {
        case LagrangianGroup::none:
            return LagrangianState::none;
        case LagrangianGroup::complete:
            return LagrangianState::complete;
        case LagrangianGroup::inInternalMesh:
            return LagrangianState::none; // <-- state is unknown
        default:
            return
                static_cast<LagrangianState>
                (
                    static_cast<label>(group)
                  - static_cast<label>(LagrangianGroup::onPatchZero)
                  + static_cast<label>(LagrangianState::onPatchZero)
                );
        case LagrangianGroup::toBeRemoved:
            return LagrangianState::toBeRemoved;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
