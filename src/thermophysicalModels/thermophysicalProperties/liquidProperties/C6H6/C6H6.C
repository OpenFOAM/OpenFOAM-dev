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

#include "C6H6.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(C6H6, 0);
    addToRunTimeSelectionTable(liquidProperties, C6H6,);
    addToRunTimeSelectionTable(liquidProperties, C6H6, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::C6H6::C6H6()
:
    liquidProperties
    (
        78.114,
        562.16,
        4.898e+6,
        0.25894,
        0.271,
        278.68,
        4.7961e+3,
        353.24,
        0.0,
        0.2108,
        1.8706e+4
    ),
    rho_("rho", 80.5511568, 0.2667, 562.16, 0.2818),
    pv_("pv", 78.05, -6275.5, -8.4443, 6.26e-06, 2),
    hl_("hl", 562.16, 649435.440510024, 0.7616, -0.5052, 0.1564, 0),
    Cp_
    (
        "Cp",
        1386.69124612745,
       -0.416058581048212,
        0.00542796425736744,
        0.0,
        0.0,
        0.0
    ),
    h_
    (
        "h",
        186141.395065592,
        1386.69124612745,
       -0.208029290524106,
        0.00180932141912248,
        0.0,
        0.0
    ),
    Cpg_
    (
        "Cpg",
        568.656066774202,
        2970.65826868423,
        1494.6,
        2203.57426325627,
        -678.15
    ),
    B_
    (
        "B",
        0.00184089919860716,
       -2.30176408838364,
       -309176.332027549,
       -5.12072099751645e+15,
       -2.90216862534245e+19
    ),
    mu_("mu", 6.764, 336.4, -2.687, 0.0, 0.0),
    mug_("mug", 3.134e-08, 0.9676, 7.9, 0.0),
    kappa_("kappa", 0.2407, -0.0003202, 0.0, 0.0, 0.0, 0.0),
    kappag_("kappag", 1.652e-05, 1.3117, 491, 0.0),
    sigma_("sigma", 562.16, 0.07195, 1.2389, 0.0, 0.0, 0.0),
    D_("D", 147.18, 20.1, 78.114, 28) // note: Same as nHeptane
{}


Foam::C6H6::C6H6
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


Foam::C6H6::C6H6(const dictionary& dict)
:
    C6H6()
{
    readIfPresent(*this, dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::C6H6::write(Ostream& os) const
{
    liquidProperties::write(*this, os);
}


// ************************************************************************* //
