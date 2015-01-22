/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2015 OpenFOAM Foundation
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

#include "fluidThermo.H"
#include "addToRunTimeSelectionTable.H"
#include "makeTurbulenceModel.H"

#include "thermalDiffusivity.H"
#include "eddyDiffusivity.H"

#include "laminar.H"
#include "RASModel.H"
#include "LESModel.H"

makeBaseTurbulenceModel
(
    geometricOneField,
    volScalarField,
    compressibleTurbulenceModel,
    thermalDiffusivity,
    fluidThermo
);

#define makeRASModel(Type)                                                     \
    makeTemplatedTurbulenceModel                                               \
    (fluidThermothermalDiffusivity, RAS, Type)

#define makeLESModel(Type)                                                     \
    makeTemplatedTurbulenceModel                                               \
    (fluidThermothermalDiffusivity, LES, Type)

#include "SpalartAllmaras.H"
makeRASModel(SpalartAllmaras);

#include "kEpsilon.H"
makeRASModel(kEpsilon);

#include "realizableKE.H"
makeRASModel(realizableKE);

#include "buoyantKEpsilon.H"
makeRASModel(buoyantKEpsilon);

#include "LaunderSharmaKE.H"
makeRASModel(LaunderSharmaKE);

#include "kOmegaSST.H"
makeRASModel(kOmegaSST);

#include "Smagorinsky.H"
makeLESModel(Smagorinsky);

#include "kEqn.H"
makeLESModel(kEqn);

#include "SpalartAllmarasDES.H"
makeLESModel(SpalartAllmarasDES);

#include "SpalartAllmarasDDES.H"
makeLESModel(SpalartAllmarasDDES);

#include "SpalartAllmarasIDDES.H"
makeLESModel(SpalartAllmarasIDDES);


// ************************************************************************* //
