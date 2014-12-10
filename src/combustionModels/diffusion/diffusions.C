/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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
#include "psiThermoCombustion.H"
#include "rhoThermoCombustion.H"
#include "diffusion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{
    // Combustion models based on sensibleEnthalpy
    makeCombustionTypesThermo
    (
        diffusion,
        psiThermoCombustion,
        gasHThermoPhysics,
        psiCombustionModel
    );

    makeCombustionTypesThermo
    (
        diffusion,
        psiThermoCombustion,
        constGasHThermoPhysics,
        psiCombustionModel
    );

    makeCombustionTypesThermo
    (
        diffusion,
        rhoThermoCombustion,
        gasHThermoPhysics,
        rhoCombustionModel
    );

    makeCombustionTypesThermo
    (
        diffusion,
        rhoThermoCombustion,
        constGasHThermoPhysics,
        rhoCombustionModel
    );

     // Combustion models based on sensibleInternalEnergy

    makeCombustionTypesThermo
    (
        diffusion,
        psiThermoCombustion,
        gasEThermoPhysics,
        psiCombustionModel
    );

    makeCombustionTypesThermo
    (
        diffusion,
        psiThermoCombustion,
        constGasEThermoPhysics,
        psiCombustionModel
    );

    makeCombustionTypesThermo
    (
        diffusion,
        rhoThermoCombustion,
        gasEThermoPhysics,
        rhoCombustionModel
    );

    makeCombustionTypesThermo
    (
        diffusion,
        rhoThermoCombustion,
        constGasEThermoPhysics,
        rhoCombustionModel
    );
}
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
