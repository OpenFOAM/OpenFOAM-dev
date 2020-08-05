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

#include "forCommonGases.H"
#include "forCommonLiquids.H"
#include "forPolynomials.H"
#include "makeChemistryReductionMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    forCommonGases(defineChemistryReductionMethod, nullArg);

    forCommonGases(makeChemistryReductionMethod, none);
    forCommonGases(makeChemistryReductionMethod, DAC);
    forCommonGases(makeChemistryReductionMethod, DRG);
    forCommonGases(makeChemistryReductionMethod, DRGEP);
    forCommonGases(makeChemistryReductionMethod, EFA);
    forCommonGases(makeChemistryReductionMethod, PFA);

    forCommonLiquids(defineChemistryReductionMethod, nullArg);

    forCommonLiquids(makeChemistryReductionMethod, none);
    forCommonLiquids(makeChemistryReductionMethod, DAC);
    forCommonLiquids(makeChemistryReductionMethod, DRG);
    forCommonLiquids(makeChemistryReductionMethod, DRGEP);
    forCommonLiquids(makeChemistryReductionMethod, EFA);
    forCommonLiquids(makeChemistryReductionMethod, PFA);

    forPolynomials(defineChemistryReductionMethod, nullArg);

    forPolynomials(makeChemistryReductionMethod, none);
    forPolynomials(makeChemistryReductionMethod, DAC);
    forPolynomials(makeChemistryReductionMethod, DRG);
    forPolynomials(makeChemistryReductionMethod, DRGEP);
    forPolynomials(makeChemistryReductionMethod, EFA);
    forPolynomials(makeChemistryReductionMethod, PFA);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
