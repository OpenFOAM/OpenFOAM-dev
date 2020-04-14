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
#include "addToRunTimeSelectionTable.H"
#include "makeThermophysicalTransportModel.H"

#include "laminarThermophysicalTransportModel.H"
#include "RASThermophysicalTransportModel.H"
#include "LESThermophysicalTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeThermophysicalTransportModelTypes
(
    PhaseThermophysicalTransportModel,
    phaseCompressibleMomentumTransportModel
);


makeThermophysicalTransportModel
(
    PhaseThermophysicalTransportModel,
    phaseCompressibleMomentumTransportModel
);


#define makeThermophysicalTransportLaminarModel(Type)                          \
    makeTemplatedThermophysicalTransportModel                                  \
    (                                                                          \
        PhaseThermophysicalTransportModel,                                     \
        phaseCompressibleMomentumTransportModel,                               \
        laminar,                                                               \
        Type                                                                   \
    )

#define makeThermophysicalTransportRASModel(Type)                              \
    makeTemplatedThermophysicalTransportModel                                  \
    (                                                                          \
        PhaseThermophysicalTransportModel,                                     \
        phaseCompressibleMomentumTransportModel,                               \
        RAS,                                                                   \
        Type                                                                   \
    )

#define makeThermophysicalTransportLESModel(Type)                              \
    makeTemplatedThermophysicalTransportModel                                  \
    (                                                                          \
        PhaseThermophysicalTransportModel,                                     \
        phaseCompressibleMomentumTransportModel,                               \
        LES,                                                                   \
        Type                                                                   \
    )


// -------------------------------------------------------------------------- //
// Laminar models
// -------------------------------------------------------------------------- //

#include "Fourier.H"
makeThermophysicalTransportLaminarModel(Fourier);


// -------------------------------------------------------------------------- //
// RAS models
// -------------------------------------------------------------------------- //

#include "RASeddyDiffusivity.H"
makeThermophysicalTransportRASModel(eddyDiffusivity);


// -------------------------------------------------------------------------- //
// LES models
// -------------------------------------------------------------------------- //

#include "LESeddyDiffusivity.H"
makeThermophysicalTransportLESModel(eddyDiffusivity);


// ************************************************************************* //
