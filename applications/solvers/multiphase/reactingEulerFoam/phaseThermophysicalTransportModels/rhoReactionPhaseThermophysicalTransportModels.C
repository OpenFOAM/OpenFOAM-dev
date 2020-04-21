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

#include "PhaseThermophysicalTransportModel.H"
#include "phaseCompressibleMomentumTransportModel.H"
#include "rhoReactionThermo.H"
#include "makeThermophysicalTransportModel.H"
#include "addToRunTimeSelectionTable.H"

#include "laminarThermophysicalTransportModel.H"
#include "RASThermophysicalTransportModel.H"
#include "LESThermophysicalTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeThermophysicalTransportModelTypes
(
    PhaseThermophysicalTransportModel,
    phaseCompressibleMomentumTransportModel,
    rhoReactionThermo
);


makeThermophysicalTransportModels
(
    PhaseThermophysicalTransportModel,
    phaseCompressibleMomentumTransportModel,
    rhoReactionThermo
);


#define makeLaminarThermophysicalTransportModel(Type)                          \
    makeThermophysicalTransportModel                                           \
    (                                                                          \
        PhaseThermophysicalTransportModel,                                     \
        phaseCompressibleMomentumTransportModel,                               \
        rhoReactionThermo,                                                     \
        laminar,                                                               \
        Type                                                                   \
    )

#define makeRASLESThermophysicalTransportModel(SType, Type)                    \
    makeTurbulenceThermophysicalTransportModel                                 \
    (                                                                          \
        PhaseThermophysicalTransportModel,                                     \
        phaseCompressibleMomentumTransportModel,                               \
        rhoReactionThermo,                                                     \
        SType,                                                                 \
        Type                                                                   \
    )

#define makeRASThermophysicalTransportModel(Type)                              \
    makeThermophysicalTransportModel                                           \
    (                                                                          \
        PhaseThermophysicalTransportModel,                                     \
        phaseCompressibleMomentumTransportModel,                               \
        rhoReactionThermo,                                                     \
        RAS,                                                                   \
        Type                                                                   \
    )

#define makeLESThermophysicalTransportModel(Type)                              \
    makeThermophysicalTransportModel                                           \
    (                                                                          \
        PhaseThermophysicalTransportModel,                                     \
        phaseCompressibleMomentumTransportModel,                               \
        rhoReactionThermo,                                                     \
        LES,                                                                   \
        Type                                                                   \
    )


// -------------------------------------------------------------------------- //
// Laminar models
// -------------------------------------------------------------------------- //

#include "Fourier.H"
makeLaminarThermophysicalTransportModel(Fourier);


// -------------------------------------------------------------------------- //
// RAS models
// -------------------------------------------------------------------------- //

#include "eddyDiffusivity.H"
makeRASLESThermophysicalTransportModel(RAS, eddyDiffusivity);

#include "nonUnityLewisEddyDiffusivity.H"
makeRASLESThermophysicalTransportModel(RAS, nonUnityLewisEddyDiffusivity);


// -------------------------------------------------------------------------- //
// LES models
// -------------------------------------------------------------------------- //

#include "eddyDiffusivity.H"
makeRASLESThermophysicalTransportModel(LES, eddyDiffusivity);

#include "nonUnityLewisEddyDiffusivity.H"
makeRASLESThermophysicalTransportModel(LES, nonUnityLewisEddyDiffusivity);


// ************************************************************************* //
