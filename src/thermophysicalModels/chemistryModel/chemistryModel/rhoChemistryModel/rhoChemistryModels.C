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

InClass
    Foam::rhoChemistryModel

Description
    Creates chemistry model instances templated on the type of thermodynamics

\*---------------------------------------------------------------------------*/

#include "makeChemistryModel.H"

#include "rhoChemistryModel.H"
#include "chemistryModel.H"
#include "thermoPhysicsTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Chemistry moldels based on sensibleEnthalpy
    makeChemistryModel
    (
        chemistryModel,
        rhoChemistryModel,
        constGasHThermoPhysics
    );

    makeChemistryModel
    (
        chemistryModel,
        rhoChemistryModel,
        gasHThermoPhysics
    );

    makeChemistryModel
    (
        chemistryModel,
        rhoChemistryModel,
        constIncompressibleGasHThermoPhysics
    );

    makeChemistryModel
    (
        chemistryModel,
        rhoChemistryModel,
        incompressibleGasHThermoPhysics
    );

    makeChemistryModel
    (
        chemistryModel,
        rhoChemistryModel,
        icoPoly8HThermoPhysics
    );


    // Chemistry moldels based on sensibleInternalEnergy
    makeChemistryModel
    (
        chemistryModel,
        rhoChemistryModel,
        constGasEThermoPhysics
    );

    makeChemistryModel
    (
        chemistryModel,
        rhoChemistryModel,
        gasEThermoPhysics
    );

    makeChemistryModel
    (
        chemistryModel,
        rhoChemistryModel,
        constIncompressibleGasEThermoPhysics
    );

    makeChemistryModel
    (
        chemistryModel,
        rhoChemistryModel,
        incompressibleGasEThermoPhysics
    );

    makeChemistryModel
    (
        chemistryModel,
        rhoChemistryModel,
        icoPoly8EThermoPhysics
    );
}

// ************************************************************************* //
