/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "C8H10.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(C8H10, 0);
    addToRunTimeSelectionTable(liquidProperties, C8H10,);
    addToRunTimeSelectionTable(liquidProperties, C8H10, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::C8H10::C8H10()
:
    liquidProperties
    (
        106.167,
        617.17,
        3.6094e+6,
        0.37381,
        0.263,
        178.15,
        4.038e-3,
        409.35,
        1.9680e-30,
        0.3036,
        1.8043e+4
    ),
    rho_("rho", 76.3765398, 0.26438, 617.17, 0.2921),
    pv_("pv", 88.246, -7691.1, -9.797, 5.931e-06, 2.0),
    hl_("hl", 617.17, 516167.924119547, 0.3882, 0.0, 0.0, 0.0),
    Cp_
    (
        "Cp",
        818.521762883005,
        6.66873887366131,
       -0.0248005500767658,
        4.23860521631015e-05,
        0.0,
        0.0
    ),
    h_
    (
        "h",
       -524002.612929508,
        818.521762883005,
        3.33436943683065,
       -0.00826685002558862,
        1.05965130407754e-05,
        0.0
    ),
    Cpg_
    (
        "Cpg",
        738.835984816374,
        3201.5598067196,
        1559,
        2285.07916772632,
        -702.0
    ),
    B_
    (
        "B",
        0.00165776559571242,
       -2.77958310962917,
       -388067.855359952,
       -5.86905535618412e+18,
        1.58052878954854e+21
    ),
    mu_("mu", -10.452, 1048.4, -0.0715, 0.0, 0.0),
    mug_("mug", 1.2e-06, 0.4518, 439.0, 0.0),
    kappa_("kappa", 0.20149, -0.00023988, 0.0, 0.0, 0.0, 0.0),
    kappag_("kappag", 1.708e-05, 1.319, 565.6, 0.0),
    sigma_("sigma", 617.17, 0.066, 1.268, 0.0, 0.0, 0.0),
    D_("D", 147.18, 20.1, 106.167, 28.0) // note: Same as nHeptane
{}


Foam::C8H10::C8H10
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
    D_(vapourDiffusivity)
{}


Foam::C8H10::C8H10(const dictionary& dict)
:
    C8H10()
{
    readIfPresent(*this, dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::C8H10::write(Ostream& os) const
{
    liquidProperties::write(*this, os);
}


// ************************************************************************* //
