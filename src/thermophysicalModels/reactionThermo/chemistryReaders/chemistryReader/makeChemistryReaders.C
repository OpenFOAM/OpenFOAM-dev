/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

makeChemistryReaderType(foamChemistryReader, constGasHThermoPhysics);
makeChemistryReaderType(foamChemistryReader, gasHThermoPhysics);
makeChemistryReaderType
(
    foamChemistryReader,
    constIncompressibleGasHThermoPhysics
);
makeChemistryReaderType(foamChemistryReader, incompressibleGasHThermoPhysics);
makeChemistryReaderType(foamChemistryReader, icoPoly8HThermoPhysics);



// Solid chemistry readers based on sensibleInternalEnergy

makeChemistryReader(constGasEThermoPhysics);
makeChemistryReader(gasEThermoPhysics);
makeChemistryReader(constIncompressibleGasEThermoPhysics);
makeChemistryReader(incompressibleGasEThermoPhysics);
makeChemistryReader(icoPoly8EThermoPhysics);

makeChemistryReaderType(foamChemistryReader, constGasEThermoPhysics);
makeChemistryReaderType(foamChemistryReader, gasEThermoPhysics);
makeChemistryReaderType
(
    foamChemistryReader,
    constIncompressibleGasEThermoPhysics
);
makeChemistryReaderType(foamChemistryReader, incompressibleGasEThermoPhysics);
makeChemistryReaderType(foamChemistryReader, icoPoly8EThermoPhysics);


// Solid chemistry readers for solids based on sensibleInternalEnergy

makeChemistryReader(hConstSolidThermoPhysics);
makeChemistryReader(hExponentialSolidThermoPhysics);
makeChemistryReader(hExpKappaConstSolidThermoPhysics);

makeChemistryReaderType(foamChemistryReader, hConstSolidThermoPhysics);
makeChemistryReaderType(foamChemistryReader, hExponentialSolidThermoPhysics);
makeChemistryReaderType(foamChemistryReader, hExpKappaConstSolidThermoPhysics);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
