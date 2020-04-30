/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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

#ifndef makeThermoCombustionModel_H
#define makeThermoCombustionModel_H

#include "CombustionModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeThermoCombustionModel(Model, ReactionThermo, ThermoPhys)           \
                                                                               \
    typedef combustionModels::Model<ReactionThermo, ThermoPhys>                \
        Model##ReactionThermo##ThermoPhys;                                     \
    typedef CombustionModel<ReactionThermo>                                    \
        Model##ReactionThermo##ThermoPhys##Base;                               \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        Model##ReactionThermo##ThermoPhys,                                     \
        (                                                                      \
            word(Model##ReactionThermo##ThermoPhys::typeName_()) + "<"         \
          + ReactionThermo::typeName + "," + ThermoPhys::typeName() + ">"      \
        ).c_str(),                                                             \
        0                                                                      \
    );                                                                         \
                                                                               \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        Model##ReactionThermo##ThermoPhys##Base,                               \
        Model##ReactionThermo##ThermoPhys,                                     \
        dictionary                                                             \
    )

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
