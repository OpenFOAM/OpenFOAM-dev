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

#include "noChemistrySolver.H"
#include "EulerImplicit.H"
#include "ode.H"

#include "StandardChemistryModel.H"
#include "TDACChemistryModel.H"

#include "psiReactionThermo.H"
#include "rhoReactionThermo.H"

#include "forCommonGases.H"
#include "forCommonLiquids.H"
#include "forPolynomials.H"
#include "makeChemistrySolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineChemistrySolvers(ReactionThermo, ThermoPhysics)                  \
    defineChemistrySolver                                                      \
    (                                                                          \
        StandardChemistryModel,                                                \
        ReactionThermo,                                                        \
        ThermoPhysics                                                          \
    );                                                                         \
    defineChemistrySolver                                                      \
    (                                                                          \
        TDACChemistryModel,                                                    \
        ReactionThermo,                                                        \
        ThermoPhysics                                                          \
    )

#define makeChemistrySolvers(Solver, ReactionThermo, ThermoPhysics)            \
    makeChemistrySolver                                                        \
    (                                                                          \
        Solver,                                                                \
        StandardChemistryModel,                                                \
        ReactionThermo,                                                        \
        ThermoPhysics                                                          \
    );                                                                         \
    makeChemistrySolver                                                        \
    (                                                                          \
        Solver,                                                                \
        TDACChemistryModel,                                                    \
        ReactionThermo,                                                        \
        ThermoPhysics                                                          \
    )

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    forCommonGases(defineChemistrySolvers, psiReactionThermo);
    forCommonGases(defineChemistrySolvers, rhoReactionThermo);

    forCommonGases(makeChemistrySolvers, noChemistrySolver, psiReactionThermo);
    forCommonGases(makeChemistrySolvers, noChemistrySolver, rhoReactionThermo);
    forCommonGases(makeChemistrySolvers, EulerImplicit, psiReactionThermo);
    forCommonGases(makeChemistrySolvers, EulerImplicit, rhoReactionThermo);
    forCommonGases(makeChemistrySolvers, ode, psiReactionThermo);
    forCommonGases(makeChemistrySolvers, ode, rhoReactionThermo);

    forCommonLiquids(defineChemistrySolvers, rhoReactionThermo);

    forCommonLiquids
    (
        makeChemistrySolvers,
        noChemistrySolver,
        rhoReactionThermo
    );
    forCommonLiquids(makeChemistrySolvers, EulerImplicit, rhoReactionThermo);
    forCommonLiquids(makeChemistrySolvers, ode, rhoReactionThermo);

    forPolynomials(defineChemistrySolvers, rhoReactionThermo);

    forPolynomials(makeChemistrySolvers, noChemistrySolver, rhoReactionThermo);
    forPolynomials(makeChemistrySolvers, EulerImplicit, rhoReactionThermo);
    forPolynomials(makeChemistrySolvers, ode, rhoReactionThermo);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
