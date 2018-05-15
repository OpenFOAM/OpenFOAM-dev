/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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

InClass
    Foam::psiChemistryModel

Description
    Creates chemistry model instances templated on the type of thermodynamics

\*---------------------------------------------------------------------------*/

#include "makeChemistryModel.H"

#include "psiReactionThermo.H"
#include "rhoReactionThermo.H"

#include "StandardChemistryModel.H"
#include "TDACChemistryModel.H"
#include "thermoPhysicsTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Make base types
    makeChemistryModel(psiReactionThermo);
    makeChemistryModel(rhoReactionThermo);

    // Chemistry moldels based on sensibleEnthalpy
    makeChemistryModelType
    (
        StandardChemistryModel,
        psiReactionThermo,
        constGasHThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        psiReactionThermo,
        gasHThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        psiReactionThermo,
        constIncompressibleGasHThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        psiReactionThermo,
        incompressibleGasHThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        psiReactionThermo,
        icoPoly8HThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        psiReactionThermo,
        constFluidHThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        psiReactionThermo,
        constAdiabaticFluidHThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        psiReactionThermo,
        constHThermoPhysics
    );


    makeChemistryModelType
    (
        StandardChemistryModel,
        rhoReactionThermo,
        constGasHThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        rhoReactionThermo,
        gasHThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        rhoReactionThermo,
        constIncompressibleGasHThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        rhoReactionThermo,
        incompressibleGasHThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        rhoReactionThermo,
        icoPoly8HThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        rhoReactionThermo,
        constFluidHThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        rhoReactionThermo,
        constAdiabaticFluidHThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        rhoReactionThermo,
        constHThermoPhysics
    );


    makeChemistryModelType
    (
        TDACChemistryModel,
        psiReactionThermo,
        constGasHThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        psiReactionThermo,
        gasHThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        psiReactionThermo,
        constIncompressibleGasHThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        psiReactionThermo,
        incompressibleGasHThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        psiReactionThermo,
        icoPoly8HThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        psiReactionThermo,
        constFluidHThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        psiReactionThermo,
        constAdiabaticFluidHThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        psiReactionThermo,
        constHThermoPhysics
    );


    makeChemistryModelType
    (
        TDACChemistryModel,
        rhoReactionThermo,
        constGasHThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        rhoReactionThermo,
        gasHThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        rhoReactionThermo,
        constIncompressibleGasHThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        rhoReactionThermo,
        incompressibleGasHThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        rhoReactionThermo,
        icoPoly8HThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        rhoReactionThermo,
        constFluidHThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        rhoReactionThermo,
        constAdiabaticFluidHThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        rhoReactionThermo,
        constHThermoPhysics
    );


    // Chemistry moldels based on sensibleInternalEnergy
    makeChemistryModelType
    (
        StandardChemistryModel,
        psiReactionThermo,
        constGasEThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        psiReactionThermo,
        gasEThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        psiReactionThermo,
        constIncompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        psiReactionThermo,
        incompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        psiReactionThermo,
        icoPoly8EThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        psiReactionThermo,
        constFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        psiReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        psiReactionThermo,
        constEThermoPhysics
    );



    makeChemistryModelType
    (
        StandardChemistryModel,
        rhoReactionThermo,
        constGasEThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        rhoReactionThermo,
        gasEThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        rhoReactionThermo,
        constIncompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        rhoReactionThermo,
        incompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        rhoReactionThermo,
        icoPoly8EThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        rhoReactionThermo,
        constFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        rhoReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        rhoReactionThermo,
        constEThermoPhysics
    );


    makeChemistryModelType
    (
        TDACChemistryModel,
        psiReactionThermo,
        constGasEThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        psiReactionThermo,
        gasEThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        psiReactionThermo,
        constIncompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        psiReactionThermo,
        incompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        psiReactionThermo,
        icoPoly8EThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        psiReactionThermo,
        constFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        psiReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        psiReactionThermo,
        constEThermoPhysics
    );


    makeChemistryModelType
    (
        TDACChemistryModel,
        rhoReactionThermo,
        constGasEThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        rhoReactionThermo,
        gasEThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        rhoReactionThermo,
        constIncompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        rhoReactionThermo,
        incompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        rhoReactionThermo,
        icoPoly8EThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        rhoReactionThermo,
        constFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        rhoReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        rhoReactionThermo,
        constEThermoPhysics
    );
}

// ************************************************************************* //
