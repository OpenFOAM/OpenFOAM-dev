/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "C16H34.H"
#include "addToRunTimeSelectionTable.H"

#include "thermodynamicConstants.H"
using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(C16H34, 0);
    addToRunTimeSelectionTable(liquidProperties, C16H34,);
    addToRunTimeSelectionTable(liquidProperties, C16H34, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::C16H34::C16H34()
:
    liquidProperties
    (
        typeName,
        226.446,
        720.60,
        1.4186e+6,
        0.93,
        0.22,
        291.32,
        8.7467e-2,
        560.01,
        0.0,
        0.7471,
        1.6052e+4
    ),
    rho_("rho", 61.94656776, 0.25442, 720.6, 0.3238),
    pv_("pv", 233.1, -17346, -32.251, 0.02407, 1.0),
    hl_("hl", 720.60, 430654.548987397, 0.4122, 0.0, 0.0, 0.0),
    Cp_
    (
        "Cp",
        3769.90540791182,
       -12.5871068599136,
        0.0247211255663602,
        0.0,
        0.0,
        0.0
    ),
    h_
    (
        "h",
       -2777201.30410301,
        3769.90540791182,
       -6.29355342995681,
        0.00824037518878673,
        0.0,
        0.0
    ),
    Cpg_
    (
        "Cpg",
        1128.74592618108,
        3600.8584828171,
        -1429.7,
        2259.69988429913,
        679.0
    ),
    B_
    (
        "B",
        0.0025091191718997,
       -2.46668079807106,
       -1704070.72767901,
       -3.00623548219001e+19,
        7.07320950690231e+21
    ),
    mu_("mu", -18.388, 2056.8, 0.98681, 0.0, 0.0),
    mug_("mug", 1.2463e-07, 0.7322, 395.0, 6000.0),
    kappa_("kappa", 0.1963, -0.00019, 0.0, 0.0, 0.0, 0.0),
    kappag_("kappag", 3.075e-06, 1.552, 678.0, 0.0),
    sigma_("sigma", 720.60, 0.05699, 1.3929, 0.0, 0.0, 0.0),
    D_("D", 147.18, 20.1, 226.446, 28.0), // note: Same as nHeptane
    hf_(h_.value(Tstd))
{}


Foam::C16H34::C16H34
(
    const liquidProperties& l,
    const Function1s::NSRDS5& density,
    const Function1s::NSRDS1& vapourPressure,
    const Function1s::NSRDS6& heatOfVapourisation,
    const Function1s::NSRDS0& heatCapacity,
    const Function1s::NSRDS0& enthalpy,
    const Function1s::NSRDS7& idealGasHeatCapacity,
    const Function1s::NSRDS4& secondVirialCoeff,
    const Function1s::NSRDS1& dynamicViscosity,
    const Function1s::NSRDS2& vapourDynamicViscosity,
    const Function1s::NSRDS0& thermalConductivity,
    const Function1s::NSRDS2& vapourThermalConductivity,
    const Function1s::NSRDS6& surfaceTension,
    const Function2s::APIdiffCoef& vapourDiffusivity
)
:
    liquidProperties(l),
    rho_(density),
    pv_(vapourPressure),
    hl_(heatOfVapourisation),
    Cp_(heatCapacity),
    h_(enthalpy),
    Cpg_(idealGasHeatCapacity),
    B_(secondVirialCoeff),
    mu_(dynamicViscosity),
    mug_(vapourDynamicViscosity),
    kappa_(thermalConductivity),
    kappag_(vapourThermalConductivity),
    sigma_(surfaceTension),
    D_(vapourDiffusivity),
    hf_(h_.value(Tstd))
{}


Foam::C16H34::C16H34(const dictionary& dict)
:
    C16H34()
{
    readIfPresent(*this, dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::C16H34::write(Ostream& os) const
{
    liquidProperties::write(*this, os);
}


// ************************************************************************* //
