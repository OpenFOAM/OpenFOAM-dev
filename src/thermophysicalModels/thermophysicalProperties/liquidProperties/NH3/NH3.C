/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
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

#include "NH3.H"
#include "thermodynamicConstants.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(NH3, 0);
    addToRunTimeSelectionTable(liquidProperties, NH3,);
    addToRunTimeSelectionTable(liquidProperties, NH3, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::NH3::NH3()
:
    liquidProperties
    (
        typeName,
        17.030,
        405.65,
        1.1278e+07,
        0.07247,
        0.242,
        195.41,
        6.1177e+03,
        239.72,
        4.9034e-30,
        0.2520,
        2.9217e+04
    ),
    rho_("rho", 3.5430e+00*W(), 2.5471e-01, Tc(), 2.8870e-01),
    pv_("pv", 9.0451e+01, -4.6690e+03, -1.1601e+01, 1.7183e-02, 1),
    hl_("hl", Tc(), 3.1523e+07/W(), 3.9140e-01, -2.2890e-01, 2.3090e-01, 0),
    Cp_
    (
        "Cp",
        1.0827e+06/W(),
        -1.5541e+04/W(),
        8.9011e+01/W(),
        -2.2513e-01/W(),
        2.1336e-04/W(),
        0
    ),
    h_(Cp_.integral("h", - Cp_.integral("", 0).value(Tstd) - 46190000/W())),
    Cpg_("Cpg", 3.3190e+04/W(), 7.4230e+04/W(), 5.4040e+02, 8.8950e-01),
    B_("B", 1.5600e-02, -1.9900e+01, -5.0500e+06, -2.5330e+18, 3.8700e+20),
    mu_("mu", -1.6430e+00, 4.5560e+02, -1.5637e+00, 0, 0),
    mug_("mug", 4.1855e-08, 9.8060e-01, 3.0800e+01, 0),
    kappa_("kappa", 1.1606e+00, -2.2840e-03, 0, 0, 0, 0),
    kappag_("kappag", -4.5900e-02, 1.6520e-01, -1.7078e+03, 0),
    sigma_("sigma", 9.1200e-02, 1.1028e+00, 0, 0, 0, 0),
    D_("D", 14.9, 20.1, W(), 28),
    hf_(h_.value(Tstd))
{}


Foam::NH3::NH3
(
    const liquidProperties& l,
    const Function1s::NSRDS5& density,
    const Function1s::NSRDS1& vapourPressure,
    const Function1s::NSRDS6& heatOfVapourisation,
    const Function1s::NSRDS0& heatCapacity,
    const Function1s::NSRDS0& enthalpy,
    const Function1s::NSRDS3& idealGasHeatCapacity,
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


Foam::NH3::NH3(const dictionary& dict)
:
    NH3()
{
    readIfPresent(*this, dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::NH3::write(Ostream& os) const
{
    liquidProperties::write(*this, os);
}


// ************************************************************************* //
