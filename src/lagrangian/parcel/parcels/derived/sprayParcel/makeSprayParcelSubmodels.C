/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "sprayCloud.H"

#include "makeParcelCloudFunctionObjects.H"

// Momentum
#include "makeParcelForces.H"
#include "makeParcelDispersionModels.H"
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

makeParcelCloudFunctionObjects(sprayCloud);

// Momentum sub-models
makeParcelForces(sprayCloud);
makeParcelDispersionModels(sprayCloud);
makeSprayParcelInjectionModels(sprayCloud);
makeParcelPatchInteractionModels(sprayCloud);
makeSprayParcelStochasticCollisionModels(sprayCloud);

// Thermo sub-models
makeParcelHeatTransferModels(sprayCloud);

// Reacting sub-models
makeReactingParcelCompositionModels(sprayCloud);
makeReactingParcelPhaseChangeModels(sprayCloud);
makeReactingParcelSurfaceFilmModels(sprayCloud);

// Spray sub-models
makeParticleForceModelType(DistortedSphereDragForce, sprayCloud);
makeSprayParcelAtomizationModels(sprayCloud);
makeSprayParcelBreakupModels(sprayCloud);


// ************************************************************************* //
