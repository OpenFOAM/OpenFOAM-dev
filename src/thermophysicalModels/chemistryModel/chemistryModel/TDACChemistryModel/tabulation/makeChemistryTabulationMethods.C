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

#include "makeChemistryTabulationMethods.H"

#include "thermoPhysicsTypes.H"

#include "psiChemistryModel.H"
#include "rhoChemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Chemistry solvers based on sensibleEnthalpy
    makeChemistryTabulationMethods(psiChemistryModel, constGasHThermoPhysics);
    makeChemistryTabulationMethods(psiChemistryModel, gasHThermoPhysics);
    makeChemistryTabulationMethods
    (
        psiChemistryModel,
        constIncompressibleGasHThermoPhysics
    );
    makeChemistryTabulationMethods
    (
        psiChemistryModel,
        incompressibleGasHThermoPhysics
    );
    makeChemistryTabulationMethods(psiChemistryModel, icoPoly8HThermoPhysics);

    makeChemistryTabulationMethods(rhoChemistryModel, constGasHThermoPhysics);

    makeChemistryTabulationMethods(rhoChemistryModel, gasHThermoPhysics);
    makeChemistryTabulationMethods
    (
        rhoChemistryModel,
        constIncompressibleGasHThermoPhysics
    );
    makeChemistryTabulationMethods
    (
        rhoChemistryModel,
        incompressibleGasHThermoPhysics
    );
    makeChemistryTabulationMethods(rhoChemistryModel, icoPoly8HThermoPhysics);

    // Chemistry solvers based on sensibleInternalEnergy

    makeChemistryTabulationMethods(psiChemistryModel, constGasEThermoPhysics);

    makeChemistryTabulationMethods(psiChemistryModel, gasEThermoPhysics);
    makeChemistryTabulationMethods
    (
        psiChemistryModel,
        constIncompressibleGasEThermoPhysics
    );
    makeChemistryTabulationMethods
    (
        psiChemistryModel,
        incompressibleGasEThermoPhysics
    );
    makeChemistryTabulationMethods(psiChemistryModel, icoPoly8EThermoPhysics);

    makeChemistryTabulationMethods(rhoChemistryModel, constGasEThermoPhysics);

    makeChemistryTabulationMethods(rhoChemistryModel, gasEThermoPhysics);
    makeChemistryTabulationMethods
    (
        rhoChemistryModel,
        constIncompressibleGasEThermoPhysics
    );
    makeChemistryTabulationMethods
    (
        rhoChemistryModel,
        incompressibleGasEThermoPhysics
    );
    makeChemistryTabulationMethods(rhoChemistryModel, icoPoly8EThermoPhysics);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
