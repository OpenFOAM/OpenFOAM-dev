/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

#include "psiChemistryModel.H"
#include "rhoChemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Chemistry solvers based on sensibleEnthalpy
    makeChemistryReductionMethods(psiChemistryModel, constGasHThermoPhysics);
    makeChemistryReductionMethods(psiChemistryModel, gasHThermoPhysics);
    makeChemistryReductionMethods
    (
        psiChemistryModel,
        constIncompressibleGasHThermoPhysics
    );
    makeChemistryReductionMethods
    (
        psiChemistryModel,
        incompressibleGasHThermoPhysics
    );
    makeChemistryReductionMethods(psiChemistryModel, icoPoly8HThermoPhysics);

    makeChemistryReductionMethods(rhoChemistryModel, constGasHThermoPhysics);
    makeChemistryReductionMethods(rhoChemistryModel, gasHThermoPhysics);
    makeChemistryReductionMethods
    (
        rhoChemistryModel,
        constIncompressibleGasHThermoPhysics
    );
    makeChemistryReductionMethods
    (
        rhoChemistryModel,
        incompressibleGasHThermoPhysics
    );
    makeChemistryReductionMethods(rhoChemistryModel, icoPoly8HThermoPhysics);


    // Chemistry solvers based on sensibleInternalEnergy
    makeChemistryReductionMethods(psiChemistryModel, constGasEThermoPhysics);
    makeChemistryReductionMethods(psiChemistryModel, gasEThermoPhysics);
    makeChemistryReductionMethods
    (
        psiChemistryModel,
        constIncompressibleGasEThermoPhysics
    );
    makeChemistryReductionMethods
    (
        psiChemistryModel,
        incompressibleGasEThermoPhysics
    );
    makeChemistryReductionMethods(psiChemistryModel, icoPoly8EThermoPhysics);

    makeChemistryReductionMethods(rhoChemistryModel, constGasEThermoPhysics);
    makeChemistryReductionMethods(rhoChemistryModel, gasEThermoPhysics);
    makeChemistryReductionMethods
    (
        rhoChemistryModel,
        constIncompressibleGasEThermoPhysics
    );
    makeChemistryReductionMethods
    (
        rhoChemistryModel,
        incompressibleGasEThermoPhysics
    );
    makeChemistryReductionMethods(rhoChemistryModel, icoPoly8EThermoPhysics);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
