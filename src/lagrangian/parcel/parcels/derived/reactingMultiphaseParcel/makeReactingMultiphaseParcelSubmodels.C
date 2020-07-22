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

#include "reactingMultiphaseCloud.H"

#include "makeParcelCloudFunctionObjects.H"

// Momentum
#include "makeParcelForces.H"
#include "makeParcelDispersionModels.H"
#include "makeReactingMultiphaseParcelInjectionModels.H" // MP variant
#include "makeParcelPatchInteractionModels.H"
#include "makeReactingMultiphaseParcelStochasticCollisionModels.H" // MP variant
#include "makeReactingParcelSurfaceFilmModels.H" // Reacting variant

// Thermodynamic
#include "makeParcelHeatTransferModels.H"

// Reacting
#include "makeReactingMultiphaseParcelCompositionModels.H" // MP Variant
#include "makeReactingParcelPhaseChangeModels.H"

// Reacting multiphase
#include "makeReactingMultiphaseParcelDevolatilisationModels.H"
#include "makeReactingMultiphaseParcelSurfaceReactionModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeParcelCloudFunctionObjects(reactingMultiphaseCloud);

// Momentum sub-models
makeParcelForces(reactingMultiphaseCloud);
makeParcelDispersionModels(reactingMultiphaseCloud);
makeReactingMultiphaseParcelInjectionModels(reactingMultiphaseCloud);
makeParcelPatchInteractionModels(reactingMultiphaseCloud);
makeReactingMultiphaseParcelStochasticCollisionModels
(
    reactingMultiphaseCloud
);
makeReactingParcelSurfaceFilmModels(reactingMultiphaseCloud);

// Thermo sub-models
makeParcelHeatTransferModels(reactingMultiphaseCloud);

// Reacting sub-models
makeReactingMultiphaseParcelCompositionModels
(
    reactingMultiphaseCloud
);
makeReactingParcelPhaseChangeModels(reactingMultiphaseCloud);

// Reacting multiphase sub-models
makeReactingMultiphaseParcelDevolatilisationModels
(
    reactingMultiphaseCloud
);
makeReactingMultiphaseParcelSurfaceReactionModels
(
    reactingMultiphaseCloud
);


// ************************************************************************* //
