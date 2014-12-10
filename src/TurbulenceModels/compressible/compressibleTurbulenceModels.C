/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "CompressibleTurbulenceModel.H"
#include "fluidThermo.H"
#include "addToRunTimeSelectionTable.H"
#include "makeTurbulenceModel.H"

#include "laminar.H"
#include "RASModel.H"
#include "LESModel.H"

makeBaseTurbulenceModel
(
    geometricOneField,
    volScalarField,
    compressibleTurbulenceModel,
    CompressibleTurbulenceModel,
    fluidThermo
);

#define makeRASModel(Type)                                                     \
    makeTemplatedTurbulenceModel                                               \
    (fluidThermoCompressibleTurbulenceModel, RAS, Type)

#define makeLESModel(Type)                                                     \
    makeTemplatedTurbulenceModel                                               \
    (fluidThermoCompressibleTurbulenceModel, LES, Type)

#include "kEpsilon.H"
makeRASModel(kEpsilon);

#include "buoyantKEpsilon.H"
makeRASModel(buoyantKEpsilon);

#include "Smagorinsky.H"
makeLESModel(Smagorinsky);

#include "kEqn.H"
makeLESModel(kEqn);


// ************************************************************************* //
