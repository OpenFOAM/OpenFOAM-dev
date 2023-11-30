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

#include "bC10H7CH3.H"
#include "addToRunTimeSelectionTable.H"

#include "thermodynamicConstants.H"
using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(bC10H7CH3, 0);
    addToRunTimeSelectionTable(liquidProperties, bC10H7CH3,);
    addToRunTimeSelectionTable(liquidProperties, bC10H7CH3, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bC10H7CH3::bC10H7CH3()
:
    liquidProperties
    (
        typeName,
        142.2,
        761.0,
        3.25e+6,
        0.507,
        0.260,
        307.73,
        1.7374e+1,
        514.20,
        1.4010e-30,
        0.3459,
        1.987e+4
    ),
    rho_("rho", 67.36014, 0.23843, 761, 0.2559),
    pv_("pv", 134.31, -12103, -16.195, 6.9659e-06, 2),
    hl_("hl", 761.0, 513150.492264416, 0.4044, 0.0, 0.0, 0.0),
    Cp_
    (
        "Cp",
        811.322081575246,
        2.30225035161744,
        0.0008628691983122,
        0.0,
        0.0,
        0.0
    ),
    h_
    (
        "h",
        45001.2311880177,
        811.322081575246,
        1.15112517580872,
        0.000287623066104079,
        0.0,
        0.0
    ),
    Cpg_
    (
        "Cpg",
        760.126582278481,
        2699.08579465542,
        1564.1,
        1994.51476793249,
        727.49
    ),
    B_
    (
        "B",
        0.00229430379746835,
       -3.53720112517581,
       -1067158.93108298,
        2.29746835443038e+18,
       -2.68438818565401e+21
    ),
    mu_("mu", -63.276, 4219, 7.5549, 0.0, 0.0),
    mug_("mug", 2.1791e-06, 0.3717, 712.53, 0.0),
    kappa_("kappa", 0.1962, -0.00018414, 0.0, 0.0, 0.0, 0.0),
    kappag_("kappag", 0.4477, -0.1282, -345.89, 2340100),
    sigma_("sigma", 761.0, 0.066442, 1.2634, 0.0, 0.0, 0.0),
    D_("D", 147.18, 20.1, 142.2, 28), // note: Same as nHeptane,
    hf_(h_.value(Tstd))
{}


Foam::bC10H7CH3::bC10H7CH3
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



Foam::bC10H7CH3::bC10H7CH3(const dictionary& dict)
:
    bC10H7CH3()
{
    readIfPresent(*this, dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::bC10H7CH3::write(Ostream& os) const
{
    liquidProperties::write(*this, os);
}


// ************************************************************************* //
