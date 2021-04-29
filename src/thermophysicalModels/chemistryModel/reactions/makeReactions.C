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
    forCoeffGases(defineReaction, nullArg);
    forCoeffLiquids(defineReaction, nullArg);


    // Irreversible/reversible/non-equilibrium-reversible reactions

    forCoeffGases(makeIRNReactions, ArrheniusReactionRate);
    forCoeffLiquids(makeIRNReactions, ArrheniusReactionRate);

    forCoeffGases(makeIRNReactions, LandauTellerReactionRate);
    forCoeffLiquids(makeIRNReactions, LandauTellerReactionRate);

    forCoeffGases(makeIRNReactions, thirdBodyArrheniusReactionRate);
    forCoeffLiquids(makeIRNReactions, thirdBodyArrheniusReactionRate);


    // Irreversible/reversible reactions

    forCoeffGases(makeIRReactions, JanevReactionRate);
    forCoeffLiquids(makeIRReactions, JanevReactionRate);

    forCoeffGases(makeIRReactions, powerSeriesReactionRate);
    forCoeffLiquids(makeIRReactions, powerSeriesReactionRate);


    // Pressure dependent reactions

    forCoeffGases
    (
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        LindemannFallOffFunction
    );
    forCoeffLiquids
    (
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        LindemannFallOffFunction
    );

    forCoeffGases
    (
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        TroeFallOffFunction
    );
    forCoeffLiquids
    (
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        TroeFallOffFunction
    );

    forCoeffGases
    (
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        SRIFallOffFunction
    );
    forCoeffLiquids
    (
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        SRIFallOffFunction
    );

    forCoeffGases
    (
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        LindemannFallOffFunction
    );
    forCoeffLiquids
    (
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        LindemannFallOffFunction
    );

    forCoeffGases
    (
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        TroeFallOffFunction
    );
    forCoeffLiquids
    (
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        TroeFallOffFunction
    );

    forCoeffGases
    (
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        SRIFallOffFunction
    );
    forCoeffLiquids
    (
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        SRIFallOffFunction
    );
}

// ************************************************************************* //
