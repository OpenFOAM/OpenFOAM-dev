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

#include "LagrangianState.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::word Foam::name(const LagrangianState state)
{
    switch (state)
    {
        case LagrangianState::none:
            return "none";
        case LagrangianState::complete:
            return "complete";
        case LagrangianState::inCell:
            return "inCell";
        case LagrangianState::onInternalFace:
            return "onInternalFace";
        default:
            return
               "onPatch"
              + name
                (
                    static_cast<label>(state)
                  - static_cast<label>(LagrangianState::onPatchZero)
                );
        case LagrangianState::toBeRemoved:
            return "toBeRemoved";
    }
}


Foam::word Foam::name(const LagrangianGroup group)
{
    switch (group)
    {
        case LagrangianGroup::none:
            return "none";
        case LagrangianGroup::complete:
            return "complete";
        case LagrangianGroup::inInternalMesh:
            return "inInternalMesh";
        default:
            return
               "onPatch"
              + name
                (
                    static_cast<label>(group)
                  - static_cast<label>(LagrangianGroup::onPatchZero)
                );
        case LagrangianGroup::toBeRemoved:
            return "toBeRemoved";
    }
}


// ************************************************************************* //
