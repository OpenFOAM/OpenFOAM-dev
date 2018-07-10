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

#include "makeReactionThermo.H"
#include "thermoPhysicsTypes.H"
#include "solidThermoPhysicsTypes.H"

#include "chemistryReader.H"
#include "foamChemistryReader.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Solid chemistry readers based on sensibleEnthalpy

makeChemistryReader(constGasHThermoPhysics);
makeChemistryReader(gasHThermoPhysics);
makeChemistryReader(constIncompressibleGasHThermoPhysics);
makeChemistryReader(incompressibleGasHThermoPhysics);
makeChemistryReader(icoPoly8HThermoPhysics);
makeChemistryReader(constFluidHThermoPhysics);
makeChemistryReader(constAdiabaticFluidHThermoPhysics);
makeChemistryReader(constHThermoPhysics);

makeChemistryReaderType(foamChemistryReader, constGasHThermoPhysics);
makeChemistryReaderType(foamChemistryReader, gasHThermoPhysics);
makeChemistryReaderType
(
    foamChemistryReader,
    constIncompressibleGasHThermoPhysics
);
makeChemistryReaderType(foamChemistryReader, incompressibleGasHThermoPhysics);
makeChemistryReaderType(foamChemistryReader, icoPoly8HThermoPhysics);
makeChemistryReaderType(foamChemistryReader, constFluidHThermoPhysics);
makeChemistryReaderType(foamChemistryReader, constAdiabaticFluidHThermoPhysics);
makeChemistryReaderType(foamChemistryReader, constHThermoPhysics);


// Solid chemistry readers based on sensibleInternalEnergy

makeChemistryReader(constGasEThermoPhysics);
makeChemistryReader(gasEThermoPhysics);
makeChemistryReader(constIncompressibleGasEThermoPhysics);
makeChemistryReader(incompressibleGasEThermoPhysics);
makeChemistryReader(icoPoly8EThermoPhysics);
makeChemistryReader(constFluidEThermoPhysics);
makeChemistryReader(constAdiabaticFluidEThermoPhysics);
makeChemistryReader(constEThermoPhysics);

makeChemistryReaderType(foamChemistryReader, constGasEThermoPhysics);
makeChemistryReaderType(foamChemistryReader, gasEThermoPhysics);
makeChemistryReaderType
(
    foamChemistryReader,
    constIncompressibleGasEThermoPhysics
);
makeChemistryReaderType(foamChemistryReader, incompressibleGasEThermoPhysics);
makeChemistryReaderType(foamChemistryReader, icoPoly8EThermoPhysics);
makeChemistryReaderType(foamChemistryReader, constFluidEThermoPhysics);
makeChemistryReaderType(foamChemistryReader, constAdiabaticFluidEThermoPhysics);
makeChemistryReaderType(foamChemistryReader, constEThermoPhysics);


// Solid chemistry readers for solids based on sensibleInternalEnergy

makeChemistryReader(hConstSolidThermoPhysics);
makeChemistryReader(hPowerSolidThermoPhysics);
makeChemistryReader(hExpKappaConstSolidThermoPhysics);

makeChemistryReaderType(foamChemistryReader, hConstSolidThermoPhysics);
makeChemistryReaderType(foamChemistryReader, hPowerSolidThermoPhysics);
makeChemistryReaderType(foamChemistryReader, hExpKappaConstSolidThermoPhysics);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
