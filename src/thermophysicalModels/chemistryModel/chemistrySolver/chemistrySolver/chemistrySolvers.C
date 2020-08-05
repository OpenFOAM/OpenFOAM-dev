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

#include "forCommonGases.H"
#include "forCommonLiquids.H"
#include "forPolynomials.H"
#include "makeChemistrySolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineChemistrySolvers(nullArg, ThermoPhysics)                         \
    defineChemistrySolver                                                      \
    (                                                                          \
        StandardChemistryModel,                                                \
        ThermoPhysics                                                          \
    );                                                                         \
    defineChemistrySolver                                                      \
    (                                                                          \
        TDACChemistryModel,                                                    \
        ThermoPhysics                                                          \
    )

#define makeChemistrySolvers(Solver, ThermoPhysics)                            \
    makeChemistrySolver                                                        \
    (                                                                          \
        Solver,                                                                \
        StandardChemistryModel,                                                \
        ThermoPhysics                                                          \
    );                                                                         \
    makeChemistrySolver                                                        \
    (                                                                          \
        Solver,                                                                \
        TDACChemistryModel,                                                    \
        ThermoPhysics                                                          \
    )

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    forCommonGases(defineChemistrySolvers, nullArg);

    forCommonGases(makeChemistrySolvers, noChemistrySolver);
    forCommonGases(makeChemistrySolvers, EulerImplicit);
    forCommonGases(makeChemistrySolvers, ode);

    forCommonLiquids(defineChemistrySolvers, nullArg);

    forCommonLiquids(makeChemistrySolvers, noChemistrySolver);
    forCommonLiquids(makeChemistrySolvers, EulerImplicit);
    forCommonLiquids(makeChemistrySolvers, ode);

    forPolynomials(defineChemistrySolvers, nullArg);

    forPolynomials(makeChemistrySolvers, noChemistrySolver);
    forPolynomials(makeChemistrySolvers, EulerImplicit);
    forPolynomials(makeChemistrySolvers, ode);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
