/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::psiChemistryModel

Description
    Creates chemistry model instances templated on the type of thermodynamics

\*---------------------------------------------------------------------------*/

#include "makeChemistryModel.H"

#include "psiChemistryModel.H"
#include "chemistryModel.H"
#include "TDACChemistryModel.H"
#include "thermoPhysicsTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Chemistry moldels based on sensibleEnthalpy
    makeChemistryModel
    (
        chemistryModel,
        psiChemistryModel,
        constGasHThermoPhysics
    );

    makeChemistryModel
    (
        chemistryModel,
        psiChemistryModel,
        gasHThermoPhysics
    );

    makeChemistryModel
    (
        chemistryModel,
        psiChemistryModel,
        constIncompressibleGasHThermoPhysics
    );

    makeChemistryModel
    (
        chemistryModel,
        psiChemistryModel,
        incompressibleGasHThermoPhysics
    );

    makeChemistryModel
    (
        chemistryModel,
        psiChemistryModel,
        icoPoly8HThermoPhysics
    );


    makeChemistryModel
    (
        TDACChemistryModel,
        psiChemistryModel,
        constGasHThermoPhysics
    );

    makeChemistryModel
    (
        TDACChemistryModel,
        psiChemistryModel,
        gasHThermoPhysics
    );

    makeChemistryModel
    (
        TDACChemistryModel,
        psiChemistryModel,
        constIncompressibleGasHThermoPhysics
    );

    makeChemistryModel
    (
        TDACChemistryModel,
        psiChemistryModel,
        incompressibleGasHThermoPhysics
    );

    makeChemistryModel
    (
        TDACChemistryModel,
        psiChemistryModel,
        icoPoly8HThermoPhysics
    );


    // Chemistry moldels based on sensibleInternalEnergy
    makeChemistryModel
    (
        chemistryModel,
        psiChemistryModel,
        constGasEThermoPhysics
    );

    makeChemistryModel
    (
        chemistryModel,
        psiChemistryModel,
        gasEThermoPhysics
    );

    makeChemistryModel
    (
        chemistryModel,
        psiChemistryModel,
        constIncompressibleGasEThermoPhysics
    );

    makeChemistryModel
    (
        chemistryModel,
        psiChemistryModel,
        incompressibleGasEThermoPhysics
    );

    makeChemistryModel
    (
        chemistryModel,
        psiChemistryModel,
        icoPoly8EThermoPhysics
    );


    makeChemistryModel
    (
        TDACChemistryModel,
        psiChemistryModel,
        constGasEThermoPhysics
    );

    makeChemistryModel
    (
        TDACChemistryModel,
        psiChemistryModel,
        gasEThermoPhysics
    );

    makeChemistryModel
    (
        TDACChemistryModel,
        psiChemistryModel,
        constIncompressibleGasEThermoPhysics
    );

    makeChemistryModel
    (
        TDACChemistryModel,
        psiChemistryModel,
        incompressibleGasEThermoPhysics
    );

    makeChemistryModel
    (
        TDACChemistryModel,
        psiChemistryModel,
        icoPoly8EThermoPhysics
    );
}

// ************************************************************************* //
