/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2021 OpenFOAM Foundation
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

#include "PhaseIncompressibleMomentumTransportModel.H"
#include "kinematicTransportModel.H"
#include "addToRunTimeSelectionTable.H"
#include "makeMomentumTransportModel.H"

#include "laminarModel.H"
#include "RASModel.H"
#include "LESModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeMomentumTransportModelTypes
(
    volScalarField,
    geometricOneField,
    incompressibleMomentumTransportModel,
    PhaseIncompressibleMomentumTransportModel,
    kinematicTransportModel
);

#define makeLaminarModel(Type)                                                 \
    makeTemplatedLaminarModel                                                  \
    (                                                                          \
        kinematicTransportModelPhaseIncompressibleMomentumTransportModel,      \
        laminar,                                                               \
        Type                                                                   \
    )

#define makeRASModel(Type)                                                     \
    makeTemplatedMomentumTransportModel                                        \
    (                                                                          \
        kinematicTransportModelPhaseIncompressibleMomentumTransportModel,      \
        RAS,                                                                   \
        Type                                                                   \
    )

#define makeLESModel(Type)                                                     \
    makeTemplatedMomentumTransportModel                                        \
    (                                                                          \
        kinematicTransportModelPhaseIncompressibleMomentumTransportModel,      \
        LES,                                                                   \
        Type                                                                   \
    )


// ************************************************************************* //
