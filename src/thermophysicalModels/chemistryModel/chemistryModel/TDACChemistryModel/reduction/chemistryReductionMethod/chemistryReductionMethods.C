/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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

#include "noChemistryReduction.H"
#include "DAC.H"
#include "DRG.H"
#include "DRGEP.H"
#include "EFA.H"
#include "PFA.H"

#include "psiReactionThermo.H"
#include "rhoReactionThermo.H"

#include "forCommonGases.H"
#include "forCommonLiquids.H"
#include "forPolynomials.H"
#include "makeChemistryReductionMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    forCommonGases(defineChemistryReductionMethod, psiReactionThermo);
    forCommonGases(defineChemistryReductionMethod, rhoReactionThermo);

    forCommonGases(makeChemistryReductionMethod, none, psiReactionThermo);
    forCommonGases(makeChemistryReductionMethod, none, rhoReactionThermo);
    forCommonGases(makeChemistryReductionMethod, DAC, psiReactionThermo);
    forCommonGases(makeChemistryReductionMethod, DAC, rhoReactionThermo);
    forCommonGases(makeChemistryReductionMethod, DRG, psiReactionThermo);
    forCommonGases(makeChemistryReductionMethod, DRG, rhoReactionThermo);
    forCommonGases(makeChemistryReductionMethod, DRGEP, psiReactionThermo);
    forCommonGases(makeChemistryReductionMethod, DRGEP, rhoReactionThermo);
    forCommonGases(makeChemistryReductionMethod, EFA, psiReactionThermo);
    forCommonGases(makeChemistryReductionMethod, EFA, rhoReactionThermo);
    forCommonGases(makeChemistryReductionMethod, PFA, psiReactionThermo);
    forCommonGases(makeChemistryReductionMethod, PFA, rhoReactionThermo);

    forCommonLiquids(defineChemistryReductionMethod, rhoReactionThermo);

    forCommonLiquids(makeChemistryReductionMethod, none, rhoReactionThermo);
    forCommonLiquids(makeChemistryReductionMethod, DAC, rhoReactionThermo);
    forCommonLiquids(makeChemistryReductionMethod, DRG, rhoReactionThermo);
    forCommonLiquids(makeChemistryReductionMethod, DRGEP, rhoReactionThermo);
    forCommonLiquids(makeChemistryReductionMethod, EFA, rhoReactionThermo);
    forCommonLiquids(makeChemistryReductionMethod, PFA, rhoReactionThermo);

    forPolynomials(defineChemistryReductionMethod, rhoReactionThermo);

    forPolynomials(makeChemistryReductionMethod, none, rhoReactionThermo);
    forPolynomials(makeChemistryReductionMethod, DAC, rhoReactionThermo);
    forPolynomials(makeChemistryReductionMethod, DRG, rhoReactionThermo);
    forPolynomials(makeChemistryReductionMethod, DRGEP, rhoReactionThermo);
    forPolynomials(makeChemistryReductionMethod, EFA, rhoReactionThermo);
    forPolynomials(makeChemistryReductionMethod, PFA, rhoReactionThermo);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
