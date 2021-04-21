/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2021 OpenFOAM Foundation
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

#include "forGases.H"
#include "forLiquids.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
const char* const Foam::Tuple2<Foam::word, Foam::scalar>::typeName
(
    "Tuple2<word,scalar>"
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    forGases(defineReaction, nullArg);
    forLiquids(defineReaction, nullArg);


    // Irreversible/reversible/non-equilibrium-reversible reactions

    forGases(makeIRNReactions, ArrheniusReactionRate);
    forLiquids(makeIRNReactions, ArrheniusReactionRate);

    forGases(makeIRNReactions, LandauTellerReactionRate);
    forLiquids(makeIRNReactions, LandauTellerReactionRate);

    forGases(makeIRNReactions, thirdBodyArrheniusReactionRate);
    forLiquids(makeIRNReactions, thirdBodyArrheniusReactionRate);


    // Irreversible/reversible reactions

    forGases(makeIRReactions, JanevReactionRate);
    forLiquids(makeIRReactions, JanevReactionRate);

    forGases(makeIRReactions, powerSeriesReactionRate);
    forLiquids(makeIRReactions, powerSeriesReactionRate);


    // Pressure dependent reactions

    forGases
    (
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        LindemannFallOffFunction
    );
    forLiquids
    (
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        LindemannFallOffFunction
    );

    forGases
    (
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        TroeFallOffFunction
    );
    forLiquids
    (
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        TroeFallOffFunction
    );

    forGases
    (
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        SRIFallOffFunction
    );
    forLiquids
    (
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        SRIFallOffFunction
    );

    forGases
    (
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        LindemannFallOffFunction
    );
    forLiquids
    (
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        LindemannFallOffFunction
    );

    forGases
    (
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        TroeFallOffFunction
    );
    forLiquids
    (
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        TroeFallOffFunction
    );

    forGases
    (
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        SRIFallOffFunction
    );
    forLiquids
    (
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        SRIFallOffFunction
    );
}

// ************************************************************************* //
