/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
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

#include "makeChemistrySolver.H"
#include "${method}ChemistryModel.H"
#include "${solver}.H"

#include "typedefThermo.H"

#include "${specie}.H"

#include "thermo.H"

// EoS
#include "${equationOfState}.H"

// Thermo
#include "${thermo}Thermo.H"
#include "${energy}.H"

// Transport
#include "${transport}Transport.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = ${SHA1sum}
    //
    // Unique function name that can be checked if the correct library version
    // has been loaded
    void ${typeName}_${SHA1sum}(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define ThermoPhysics                                                          \
    ${transport}Transport${energy}${thermo}Thermo${equationOfState}${specie}

namespace Foam
{
    typedefThermo
    (
        ${transport}Transport,
        ${energy},
        ${thermo}Thermo,
        ${equationOfState},
        ${specie}
    );

    defineChemistrySolver(${method}ChemistryModel, ThermoPhysics);
    makeChemistrySolver(${solver}, ${method}ChemistryModel, ThermoPhysics);
}

#define standardChemistryModelCppTest 0
#define TDACChemistryModelCppTest 1

#if ${method}CppTest == TDACChemistryModelCppTest

#include "makeChemistryReductionMethod.H"
#include "makeChemistryTabulationMethod.H"

namespace Foam
{
    defineChemistrySolver(standardChemistryModel, ThermoPhysics);
    defineChemistryReductionMethod(nullArg, ThermoPhysics);
    defineChemistryTabulationMethod(nullArg, ThermoPhysics);
}

#include "noChemistryReduction.H"
#include "DAC.H"
#include "DRG.H"
#include "DRGEP.H"
#include "EFA.H"
#include "PFA.H"

namespace Foam
{
    makeChemistryReductionMethod(none, ThermoPhysics);
    makeChemistryReductionMethod(DAC, ThermoPhysics);
    makeChemistryReductionMethod(DRG, ThermoPhysics);
    makeChemistryReductionMethod(DRGEP, ThermoPhysics);
    makeChemistryReductionMethod(EFA, ThermoPhysics);
    makeChemistryReductionMethod(PFA, ThermoPhysics);
}

#include "noChemistryTabulation.H"
#include "ISAT.H"
namespace Foam
{
    makeChemistryTabulationMethod(none, ThermoPhysics);
    makeChemistryTabulationMethod(ISAT, ThermoPhysics);
}

#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "makeReaction.H"

#include "ArrheniusReactionRate.H"
#include "LandauTellerReactionRate.H"
#include "thirdBodyArrheniusReactionRate.H"

#include "JanevReactionRate.H"
#include "powerSeriesReactionRate.H"

#include "ChemicallyActivatedReactionRate.H"
#include "FallOffReactionRate.H"

#include "LindemannFallOffFunction.H"
#include "SRIFallOffFunction.H"
#include "TroeFallOffFunction.H"

#include "LangmuirHinshelwoodReactionRate.H"

#include "MichaelisMentenReactionRate.H"

namespace Foam
{
    defineReaction(nullArg, ThermoPhysics);

    makeIRNReactions(ArrheniusReactionRate, ThermoPhysics);
    makeIRNReactions(LandauTellerReactionRate, ThermoPhysics);
    makeIRNReactions(thirdBodyArrheniusReactionRate, ThermoPhysics);
    makeIRReactions(JanevReactionRate, ThermoPhysics);
    makeIRReactions(powerSeriesReactionRate, ThermoPhysics);

    makeIRRPressureDependentReactions
    (
        FallOffReactionRate,
        ArrheniusReactionRate,
        LindemannFallOffFunction,
        ThermoPhysics
    );

    makeIRRPressureDependentReactions
    (
        FallOffReactionRate,
        ArrheniusReactionRate,
        TroeFallOffFunction,
        ThermoPhysics
    );

    makeIRRPressureDependentReactions
    (
        FallOffReactionRate,
        ArrheniusReactionRate,
        SRIFallOffFunction,
        ThermoPhysics
    );

    makeIRRPressureDependentReactions
    (
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        LindemannFallOffFunction,
        ThermoPhysics
    );

    makeIRRPressureDependentReactions
    (
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        TroeFallOffFunction,
        ThermoPhysics
    );

    makeIRRPressureDependentReactions
    (
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        SRIFallOffFunction,
        ThermoPhysics
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "LangmuirHinshelwoodReactionRate.H"

namespace Foam
{
    makeIRReactions(LangmuirHinshelwoodReactionRate, ThermoPhysics);
}


#include "fluxLimitedLangmuirHinshelwoodReactionRate.H"

namespace Foam
{
    makeGeneralReaction
    (
        IrreversibleReaction,
        fluxLimitedLangmuirHinshelwoodReactionRate,
        ThermoPhysics
    );
}


// ************************************************************************* //
