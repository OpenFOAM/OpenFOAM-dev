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

#include "makeCombustionTypes.H"

#include "thermoPhysicsTypes.H"
#include "psiCombustionModel.H"
#include "psiThermoCombustion.H"
#include "rhoCombustionModel.H"
#include "rhoThermoCombustion.H"
#include "FSD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{
    // Combustion models based on sensibleEnthalpy
    makeCombustionTypesThermo
    (
        FSD,
        psiThermoCombustion,
        gasHThermoPhysics,
        psiCombustionModel
    );

    makeCombustionTypesThermo
    (
        FSD,
        psiThermoCombustion,
        constGasHThermoPhysics,
        psiCombustionModel
    );

    makeCombustionTypesThermo
    (
        FSD,
        rhoThermoCombustion,
        gasHThermoPhysics,
        rhoCombustionModel
    );

    makeCombustionTypesThermo
    (
        FSD,
        rhoThermoCombustion,
        constGasHThermoPhysics,
        rhoCombustionModel
    );

    // Combustion models based on sensibleInternalEnergy
    makeCombustionTypesThermo
    (
        FSD,
        psiThermoCombustion,
        gasEThermoPhysics,
        psiCombustionModel
    );

    makeCombustionTypesThermo
    (
        FSD,
        psiThermoCombustion,
        constGasEThermoPhysics,
        psiCombustionModel
    );

    makeCombustionTypesThermo
    (
        FSD,
        rhoThermoCombustion,
        gasEThermoPhysics,
        rhoCombustionModel
    );

    makeCombustionTypesThermo
    (
        FSD,
        rhoThermoCombustion,
        constGasEThermoPhysics,
        rhoCombustionModel
    );
}
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
