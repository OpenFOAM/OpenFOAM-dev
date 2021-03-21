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

#include "iC3H8O.H"
#include "addToRunTimeSelectionTable.H"

#include "thermodynamicConstants.H"
using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(iC3H8O, 0);
    addToRunTimeSelectionTable(liquidProperties, iC3H8O,);
    addToRunTimeSelectionTable(liquidProperties, iC3H8O, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::iC3H8O::iC3H8O()
:
    liquidProperties
    (
        typeName,
        60.096,
        508.31,
        4.7643e+6,
        0.22013,
        0.248,
        185.28,
        3.20e-2,
        355.41,
        5.5372e-30,
        0.6689,
        2.3575e+4
    ),
    rho_("rho", 70.91328, 0.26475, 508.31, 0.243),
    pv_("pv", 92.935, -8177.1, -10.031, 3.9988e-06, 2.0),
    hl_("hl", 508.31, 948149.627263046, 0.087, 0.3007, 0.0, 0.0),
    Cp_
    (
        "Cp",
        7760.91586794462,
       -68.3672790202343,
        0.241380457933972,
       -0.000235057241746539,
        0.0,
        0.0
    ),
    h_
    (
        "h",
       -6227786.27583977,
        7760.91586794462,
       -34.1836395101172,
        0.0804601526446574,
       -5.87643104366347e-05,
        0.0
    ),
    Cpg_
    (
        "Cpg",
        789.73642172524,
        3219.8482428115,
        1124,
        1560.83599574015,
        460.0
    ),
    B_
    (
        "B",
        0.000502529286474973,
       -0.104665867944622,
       -717185.83599574,
        3.3047124600639e+18,
       -1.43270766773163e+21
    ),
    mu_("mu", -8.23, 2282.2, -0.98495, 0.0, 0.0),
    mug_("mug", 1.993e-07, 0.7233, 178.0, 0.0),
    kappa_("kappa", 0.2029, -0.0002278, 0.0, 0.0, 0.0, 0.0),
    kappag_("kappag", -80.642, -1.4549, -604.42, 0.0),
    sigma_("sigma", 0.03818, -3.818e-05, -6.51e-08, 0.0, 0.0, 0.0),
    D_("D", 4.75e-10, 1.75, 0.0, 0.0, 0.0),
    Hf_(h_.value(Tstd))
{}


Foam::iC3H8O::iC3H8O
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
    const Function1s::NSRDS0& surfaceTension,
    const Function1s::NSRDS1& vapourDiffusivity
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


Foam::iC3H8O::iC3H8O(const dictionary& dict)
:
    iC3H8O()
{
    readIfPresent(*this, dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::iC3H8O::write(Ostream& os) const
{
    liquidProperties::write(*this, os);
}


// ************************************************************************* //
