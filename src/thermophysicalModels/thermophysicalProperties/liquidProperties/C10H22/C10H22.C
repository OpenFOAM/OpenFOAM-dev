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

#include "C10H22.H"
#include "addToRunTimeSelectionTable.H"

#include "thermodynamicConstants.H"
using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(C10H22, 0);
    addToRunTimeSelectionTable(liquidProperties, C10H22,);
    addToRunTimeSelectionTable(liquidProperties, C10H22, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::C10H22::C10H22()
:
    liquidProperties
    (
        typeName,
        142.285,
        617.70,
        2.11e+6,
        0.6,
        0.247,
        243.51,
        1.393,
        447.30,
        0.0,
        0.4923,
        1.57e+4
    ),
    rho_("rho", 60.94208835, 0.25745, 617.7, 0.28912),
    pv_("pv", 112.73, -9749.6, -13.245, 7.1266e-06, 2.0),
    hl_("hl", 617.70, 464743.296904101, 0.39797, 0.0, 0.0, 0.0),
    Cp_
    (
        "Cp",
        1958.18252099659,
       -1.39094071757388,
        0.00754612221948905,
        0.0,
        0.0,
        0.0
    ),
    h_
    (
        "h",
       -2699436.15229142,
        1958.18252099659,
       -0.695470358786942,
        0.00251537407316302,
        0.0,
        0.0
    ),
    Cpg_
    (
        "Cpg",
        1175.10630073444,
        3762.16748076045,
        1614.1,
        2658.04547211582,
        742
    ),
    B_
    (
        "B",
        0.00337351091119935,
       -4.13606494008504,
       -534560.916470464,
       -1.13364022911762e+19,
        2.80704220402713e+21
    ),
    mu_("mu", -16.468, 1533.5, 0.7511, 0.0, 0.0),
    mug_("mug", 2.64e-08, 0.9487, 71.0, 0.0),
    kappa_("kappa", 0.2063, -0.000254, 0.0, 0.0, 0.0, 0.0),
    kappag_("kappag", -668.4, 0.9323, -4071000000.0, 0.0),
    sigma_("sigma", 617.70, 0.055435, 1.3095, 0.0, 0.0, 0.0),
    D_("D", 147.18, 20.1, 142.285, 28.0), // note: Same as nHeptane
    Hf_(h_.value(Tstd))
{}


Foam::C10H22::C10H22
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


Foam::C10H22::C10H22(const dictionary& dict)
:
    C10H22()
{
    readIfPresent(*this, dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::C10H22::write(Ostream& os) const
{
    liquidProperties::write(*this, os);
}


// ************************************************************************* //
