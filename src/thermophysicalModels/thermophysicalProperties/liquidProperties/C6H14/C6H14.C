/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "C6H14.H"
#include "addToRunTimeSelectionTable.H"

#include "thermodynamicConstants.H"
using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(C6H14, 0);
    addToRunTimeSelectionTable(liquidProperties, C6H14,);
    addToRunTimeSelectionTable(liquidProperties, C6H14, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::C6H14::C6H14()
:
    liquidProperties
    (
        typeName,
        86.177,
        507.60,
        3.025e+6,
        0.371,
        0.266,
        177.83,
        9.017e-1,
        341.88,
        0.0,
        0.3013,
        1.49e+4
    ),
    rho_("rho", 61.03399848, 0.26411, 507.6, 0.27537),
    pv_("pv", 104.65, -6995.5, -12.702, 1.2381e-05, 2.0),
    hl_("hl", 507.60, 527286.863084118, 0.39002, 0.0, 0.0, 0.0),
    Cp_
    (
        "Cp",
        1997.28465831951,
       -2.13258758137322,
        0.0102964828202421,
        0.0,
        0.0,
        0.0
    ),
    h_
    (
        "h",
       -2902186.5403246,
        1997.28465831951,
       -1.06629379068661,
        0.00343216094008069,
        0.0,
        0.0
    ),
    Cpg_
    (
        "Cpg",
        1211.4601343746,
        4088.0977523005,
        1694.6,
        2748.99335089409,
        761.6
    ),
    B_
    (
        "B",
        0.0022859927822969,
       -2.32080485512376,
       -430509.300625457,
        1.93787205402834e+17,
       -7.17128700233241e+19
    ),
    mu_("mu", -20.715, 1207.5, 1.4993, 0.0, 0.0),
    mug_("mug", 1.7514e-07, 0.70737, 157.14, 0.0),
    kappa_("kappa", 0.22492, -0.0003533, 0.0, 0.0, 0.0, 0.0),
    kappag_("kappag", -650.5, 0.8053, -1412100000, 0.0),
    sigma_("sigma", 507.60, 0.055003, 1.2674, 0.0, 0.0, 0.0),
    D_("D", 147.18, 20.1, 86.177, 28), // note: Same as nHeptane
    Hf_(h_.value(Tstd))
{}


Foam::C6H14::C6H14
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
    Hf_(h_.value(Tstd))
{}


Foam::C6H14::C6H14(const dictionary& dict)
:
    C6H14()
{
    readIfPresent(*this, dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::C6H14::write(Ostream& os) const
{
    liquidProperties::write(*this, os);
}


// ************************************************************************* //
