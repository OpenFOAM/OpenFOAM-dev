/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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

#include "makeChemistryReductionMethods.H"

#include "thermoPhysicsTypes.H"

#include "psiReactionThermo.H"
#include "rhoReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Chemistry solvers based on sensibleEnthalpy
    makeChemistryReductionMethods(psiReactionThermo, constGasHThermoPhysics);
    makeChemistryReductionMethods(psiReactionThermo, gasHThermoPhysics);
    makeChemistryReductionMethods
    (
        psiReactionThermo,
        constIncompressibleGasHThermoPhysics
    );
    makeChemistryReductionMethods
    (
        psiReactionThermo,
        incompressibleGasHThermoPhysics
    );
    makeChemistryReductionMethods(psiReactionThermo, icoPoly8HThermoPhysics);
    makeChemistryReductionMethods(psiReactionThermo, constFluidHThermoPhysics);
    makeChemistryReductionMethods
    (
        psiReactionThermo,
        constAdiabaticFluidHThermoPhysics
    );
    makeChemistryReductionMethods(psiReactionThermo, constHThermoPhysics);

    makeChemistryReductionMethods(rhoReactionThermo, constGasHThermoPhysics);
    makeChemistryReductionMethods(rhoReactionThermo, gasHThermoPhysics);
    makeChemistryReductionMethods
    (
        rhoReactionThermo,
        constIncompressibleGasHThermoPhysics
    );
    makeChemistryReductionMethods
    (
        rhoReactionThermo,
        incompressibleGasHThermoPhysics
    );
    makeChemistryReductionMethods(rhoReactionThermo, icoPoly8HThermoPhysics);
    makeChemistryReductionMethods(rhoReactionThermo, constFluidHThermoPhysics);
    makeChemistryReductionMethods
    (
        rhoReactionThermo,
        constAdiabaticFluidHThermoPhysics
    );
    makeChemistryReductionMethods(rhoReactionThermo, constHThermoPhysics);


    // Chemistry solvers based on sensibleInternalEnergy
    makeChemistryReductionMethods(psiReactionThermo, constGasEThermoPhysics);
    makeChemistryReductionMethods(psiReactionThermo, gasEThermoPhysics);
    makeChemistryReductionMethods
    (
        psiReactionThermo,
        constIncompressibleGasEThermoPhysics
    );
    makeChemistryReductionMethods
    (
        psiReactionThermo,
        incompressibleGasEThermoPhysics
    );
    makeChemistryReductionMethods(psiReactionThermo, icoPoly8EThermoPhysics);
    makeChemistryReductionMethods(psiReactionThermo, constFluidEThermoPhysics);
    makeChemistryReductionMethods
    (
        psiReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );
    makeChemistryReductionMethods(psiReactionThermo, constEThermoPhysics);

    makeChemistryReductionMethods(rhoReactionThermo, constGasEThermoPhysics);
    makeChemistryReductionMethods(rhoReactionThermo, gasEThermoPhysics);
    makeChemistryReductionMethods
    (
        rhoReactionThermo,
        constIncompressibleGasEThermoPhysics
    );
    makeChemistryReductionMethods
    (
        rhoReactionThermo,
        incompressibleGasEThermoPhysics
    );
    makeChemistryReductionMethods(rhoReactionThermo, icoPoly8EThermoPhysics);
    makeChemistryReductionMethods(rhoReactionThermo, constFluidEThermoPhysics);
    makeChemistryReductionMethods
    (
        rhoReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );
    makeChemistryReductionMethods(rhoReactionThermo, constEThermoPhysics);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
