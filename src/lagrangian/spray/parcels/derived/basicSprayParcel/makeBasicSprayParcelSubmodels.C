/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "basicSprayCloud.H"

#include "makeParcelCloudFunctionObjects.H"

// Kinematic
#include "makeThermoParcelForces.H" // thermo variant
#include "makeThermoParcelTurbulenceForces.H" // add turbulence variant
#include "makeParcelDispersionModels.H"
#include "makeParcelTurbulenceDispersionModels.H" // add turbulence variant
#include "makeSprayParcelInjectionModels.H" // Spray variant
#include "makeParcelPatchInteractionModels.H"
#include "makeSprayParcelStochasticCollisionModels.H" // Spray variant

// Thermodynamic
#include "makeParcelHeatTransferModels.H"

// Reacting
#include "makeReactingParcelCompositionModels.H"
#include "makeReactingParcelPhaseChangeModels.H"
#include "makeReactingParcelSurfaceFilmModels.H"

// Spray
#include "DistortedSphereDragForce.H"
#include "makeSprayParcelAtomizationModels.H"
#include "makeSprayParcelBreakupModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeParcelCloudFunctionObjects(basicSprayCloud);

// Kinematic sub-models
makeThermoParcelForces(basicSprayCloud);
makeThermoParcelTurbulenceForces(basicSprayCloud);
makeParcelDispersionModels(basicSprayCloud);
makeParcelTurbulenceDispersionModels(basicSprayCloud);
makeSprayParcelInjectionModels(basicSprayCloud);
makeParcelPatchInteractionModels(basicSprayCloud);
makeSprayParcelStochasticCollisionModels(basicSprayCloud);

// Thermo sub-models
makeParcelHeatTransferModels(basicSprayCloud);

// Reacting sub-models
makeReactingParcelCompositionModels(basicSprayCloud);
makeReactingParcelPhaseChangeModels(basicSprayCloud);
makeReactingParcelSurfaceFilmModels(basicSprayCloud);

// Spray sub-models
makeParticleForceModelType(DistortedSphereDragForce, basicSprayCloud);
makeSprayParcelAtomizationModels(basicSprayCloud);
makeSprayParcelBreakupModels(basicSprayCloud);


// ************************************************************************* //
