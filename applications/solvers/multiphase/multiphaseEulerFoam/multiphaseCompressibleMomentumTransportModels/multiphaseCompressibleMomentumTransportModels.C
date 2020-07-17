/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2020 OpenFOAM Foundation
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

#include "phaseCompressibleMomentumTransportModel.H"
#include "addToRunTimeSelectionTable.H"
#include "makeMomentumTransportModel.H"

#include "laminarModel.H"
#include "RASModel.H"
#include "LESModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeMomentumTransportModelTypes
(
    volScalarField,
    volScalarField,
    compressibleMomentumTransportModel,
    PhaseCompressibleMomentumTransportModel,
    phaseModel
);

makeBaseMomentumTransportModel
(
    volScalarField,
    volScalarField,
    compressibleMomentumTransportModel,
    PhaseCompressibleMomentumTransportModel,
    phaseModel
);

#define makeLaminarModel(Type)                                                 \
    makeTemplatedLaminarModel                                                  \
    (phaseModelPhaseCompressibleMomentumTransportModel, laminar, Type)

#define makeRASModel(Type)                                                     \
    makeTemplatedMomentumTransportModel                                        \
    (phaseModelPhaseCompressibleMomentumTransportModel, RAS, Type)

#define makeLESModel(Type)                                                     \
    makeTemplatedMomentumTransportModel                                        \
    (phaseModelPhaseCompressibleMomentumTransportModel, LES, Type)

#include "Stokes.H"
makeLaminarModel(Stokes);

#include "generalizedNewtonian.H"
makeLaminarModel(generalizedNewtonian);

#include "kEpsilon.H"
makeRASModel(kEpsilon);

#include "LaheyKEpsilon.H"
makeRASModel(LaheyKEpsilon);

#include "kOmegaSST.H"
makeRASModel(kOmegaSST);

#include "kOmegaSSTSato.H"
makeRASModel(kOmegaSSTSato);

#include "continuousGasKEpsilon.H"
makeRASModel(continuousGasKEpsilon);

#include "mixtureKEpsilon.H"
makeRASModel(mixtureKEpsilon);

#include "Smagorinsky.H"
makeLESModel(Smagorinsky);

#include "kEqn.H"
makeLESModel(kEqn);

#include "SmagorinskyZhang.H"
makeLESModel(SmagorinskyZhang);

#include "NicenoKEqn.H"
makeLESModel(NicenoKEqn);

#include "continuousGasKEqn.H"
makeLESModel(continuousGasKEqn);

#include "kineticTheoryModel.H"
makeMomentumTransportModel
(phaseModelPhaseCompressibleMomentumTransportModel, RAS, kineticTheoryModel);

#include "phasePressureModel.H"
makeMomentumTransportModel
(phaseModelPhaseCompressibleMomentumTransportModel, RAS, phasePressureModel);


// ************************************************************************* //
