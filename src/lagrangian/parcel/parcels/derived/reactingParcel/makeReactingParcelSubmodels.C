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

#include "reactingCloud.H"

#include "makeParcelCloudFunctionObjects.H"

// Momentum
#include "makeParcelForces.H"
#include "makeParcelDispersionModels.H"
#include "makeReactingParcelInjectionModels.H" // Reacting variant
#include "makeParcelPatchInteractionModels.H"
#include "makeParcelStochasticCollisionModels.H"
#include "makeReactingParcelSurfaceFilmModels.H" // Reacting variant

// Thermodynamic
#include "makeParcelHeatTransferModels.H"

// Reacting
#include "makeReactingParcelCompositionModels.H"
#include "makeReactingParcelPhaseChangeModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeParcelCloudFunctionObjects(reactingCloud);

// Momentum sub-models
makeParcelForces(reactingCloud);
makeParcelDispersionModels(reactingCloud);
makeReactingParcelInjectionModels(reactingCloud);
makeParcelPatchInteractionModels(reactingCloud);
makeParcelStochasticCollisionModels(reactingCloud);
makeReactingParcelSurfaceFilmModels(reactingCloud);

// Thermo sub-models
makeParcelHeatTransferModels(reactingCloud);

// Reacting sub-models
makeReactingParcelCompositionModels(reactingCloud);
makeReactingParcelPhaseChangeModels(reactingCloud);


// ************************************************************************* //
