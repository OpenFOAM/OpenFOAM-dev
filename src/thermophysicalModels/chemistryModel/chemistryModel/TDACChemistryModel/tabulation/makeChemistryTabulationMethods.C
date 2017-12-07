/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenFOAM Foundation
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

#include "makeChemistryTabulationMethods.H"

#include "thermoPhysicsTypes.H"

#include "psiReactionThermo.H"
#include "rhoReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Chemistry solvers based on sensibleEnthalpy
    makeChemistryTabulationMethods(psiReactionThermo, constGasHThermoPhysics);
    makeChemistryTabulationMethods(psiReactionThermo, gasHThermoPhysics);
    makeChemistryTabulationMethods
    (
        psiReactionThermo,
        constIncompressibleGasHThermoPhysics
    );
    makeChemistryTabulationMethods
    (
        psiReactionThermo,
        incompressibleGasHThermoPhysics
    );
    makeChemistryTabulationMethods(psiReactionThermo, icoPoly8HThermoPhysics);

    makeChemistryTabulationMethods(rhoReactionThermo, constGasHThermoPhysics);

    makeChemistryTabulationMethods(rhoReactionThermo, gasHThermoPhysics);
    makeChemistryTabulationMethods
    (
        rhoReactionThermo,
        constIncompressibleGasHThermoPhysics
    );
    makeChemistryTabulationMethods
    (
        rhoReactionThermo,
        incompressibleGasHThermoPhysics
    );
    makeChemistryTabulationMethods(rhoReactionThermo, icoPoly8HThermoPhysics);

    // Chemistry solvers based on sensibleInternalEnergy

    makeChemistryTabulationMethods(psiReactionThermo, constGasEThermoPhysics);

    makeChemistryTabulationMethods(psiReactionThermo, gasEThermoPhysics);
    makeChemistryTabulationMethods
    (
        psiReactionThermo,
        constIncompressibleGasEThermoPhysics
    );
    makeChemistryTabulationMethods
    (
        psiReactionThermo,
        incompressibleGasEThermoPhysics
    );
    makeChemistryTabulationMethods(psiReactionThermo, icoPoly8EThermoPhysics);

    makeChemistryTabulationMethods(rhoReactionThermo, constGasEThermoPhysics);

    makeChemistryTabulationMethods(rhoReactionThermo, gasEThermoPhysics);
    makeChemistryTabulationMethods
    (
        rhoReactionThermo,
        constIncompressibleGasEThermoPhysics
    );
    makeChemistryTabulationMethods
    (
        rhoReactionThermo,
        incompressibleGasEThermoPhysics
    );
    makeChemistryTabulationMethods(rhoReactionThermo, icoPoly8EThermoPhysics);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
