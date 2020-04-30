/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2020 OpenFOAM Foundation
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
#include "infiniteReactionRate.H"
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

#include "forCommonGases.H"
#include "forCommonLiquids.H"
#include "forPolynomials.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
const char* const Foam::Tuple2<Foam::word, Foam::scalar>::typeName
(
    "Tuple2<word,scalar>"
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    forCommonGases(defineReaction, nullArg);
    forCommonLiquids(defineReaction, nullArg);
    forPolynomials(defineReaction, nullArg);


    // Irreversible/reversible/non-equilibrium-reversible reactions

    forCommonGases(makeIRNReactions, ArrheniusReactionRate);
    forCommonLiquids(makeIRNReactions, ArrheniusReactionRate);
    forPolynomials(makeIRNReactions, ArrheniusReactionRate);

    forCommonGases(makeIRNReactions, infiniteReactionRate);
    forCommonLiquids(makeIRNReactions, infiniteReactionRate);
    forPolynomials(makeIRNReactions, infiniteReactionRate);

    forCommonGases(makeIRNReactions, LandauTellerReactionRate);
    forCommonLiquids(makeIRNReactions, LandauTellerReactionRate);
    forPolynomials(makeIRNReactions, LandauTellerReactionRate);

    forCommonGases(makeIRNReactions, thirdBodyArrheniusReactionRate);
    forCommonLiquids(makeIRNReactions, thirdBodyArrheniusReactionRate);
    forPolynomials(makeIRNReactions, thirdBodyArrheniusReactionRate);


    // Irreversible/reversible reactions

    forCommonGases(makeIRReactions, JanevReactionRate);
    forCommonLiquids(makeIRReactions, JanevReactionRate);
    forPolynomials(makeIRReactions, JanevReactionRate);

    forCommonGases(makeIRReactions, powerSeriesReactionRate);
    forCommonLiquids(makeIRReactions, powerSeriesReactionRate);
    forPolynomials(makeIRReactions, powerSeriesReactionRate);


    // Pressure dependent reactions

    forCommonGases
    (
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        LindemannFallOffFunction
    );
    forCommonLiquids
    (
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        LindemannFallOffFunction
    );
    forPolynomials
    (
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        LindemannFallOffFunction
    );

    forCommonGases
    (
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        TroeFallOffFunction
    );
    forCommonLiquids
    (
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        TroeFallOffFunction
    );
    forPolynomials
    (
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        TroeFallOffFunction
    );

    forCommonGases
    (
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        SRIFallOffFunction
    );
    forCommonLiquids
    (
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        SRIFallOffFunction
    );
    forPolynomials
    (
        makeIRRPressureDependentReactions,
        FallOffReactionRate,
        ArrheniusReactionRate,
        SRIFallOffFunction
    );

    forCommonGases
    (
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        LindemannFallOffFunction
    );
    forCommonLiquids
    (
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        LindemannFallOffFunction
    );
    forPolynomials
    (
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        LindemannFallOffFunction
    );

    forCommonGases
    (
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        TroeFallOffFunction
    );
    forCommonLiquids
    (
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        TroeFallOffFunction
    );
    forPolynomials
    (
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        TroeFallOffFunction
    );

    forCommonGases
    (
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        SRIFallOffFunction
    );
    forCommonLiquids
    (
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        SRIFallOffFunction
    );
    forPolynomials
    (
        makeIRRPressureDependentReactions,
        ChemicallyActivatedReactionRate,
        ArrheniusReactionRate,
        SRIFallOffFunction
    );
}

// ************************************************************************* //
