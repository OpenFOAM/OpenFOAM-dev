/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "C2H5OH.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(C2H5OH, 0);
    addToRunTimeSelectionTable(liquidProperties, C2H5OH,);
    addToRunTimeSelectionTable(liquidProperties, C2H5OH, Istream);
    addToRunTimeSelectionTable(liquidProperties, C2H5OH, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::C2H5OH::C2H5OH()
:
    liquidProperties
    (
        46.069,
        516.25,
        6.3835e+6,
        0.16692,
        0.248,
        159.05,
        7.1775e-5,
        351.44,
        5.6372e-30,
        0.6371,
        2.6421e+4
    ),
    rho_(70.1308387, 0.26395, 516.25, 0.2367),
    pv_(59.796, -6595, -5.0474, 6.3e-07, 2),
    hl_(516.25, 958345.091059064, -0.4134, 0.75362, 0.0, 0.0),
    Cp_
    (
        2052.57331394213,
       -1.21990926653498,
        0.00714146172046278,
        5.20523562482363e-05,
        0.0,
        0.0
    ),
    h_
    (
       -6752827.25039109,
        2052.57331394213,
       -0.60995463326749,
        0.00238048724015426,
        1.30130890620591e-05,
        0.0
    ),
    Cpg_(909.505307256507, 3358.00646855803, 1530, 2029.56434912848, 640),
    B_
    (
       -0.00358158414552085,
        3.90718270420456,
       -1180837.43949293,
        9.81136990166923e+18,
       -3.58592545963663e+21
    ),
    mu_(8.049, 776, -3.068, 0.0, 0.0),
    mug_(1.0613e-07, 0.8066, 52.7, 0.0),
    K_(0.253, -0.000281, 0.0, 0.0, 0.0, 0.0),
    Kg_(-3.12, 0.7152, -3550000.0, 0.0),
    sigma_(3.7640e-02, -2.1570e-05, -1.025e-07, 0.0, 0.0, 0.0),
    D_(147.18, 20.1, 46.069, 28) // note: Same as nHeptane
{}


Foam::C2H5OH::C2H5OH
(
    const liquidProperties& l,
    const NSRDSfunc5& density,
    const NSRDSfunc1& vapourPressure,
    const NSRDSfunc6& heatOfVapourisation,
    const NSRDSfunc0& heatCapacity,
    const NSRDSfunc0& enthalpy,
    const NSRDSfunc7& idealGasHeatCapacity,
    const NSRDSfunc4& secondVirialCoeff,
    const NSRDSfunc1& dynamicViscosity,
    const NSRDSfunc2& vapourDynamicViscosity,
    const NSRDSfunc0& thermalConductivity,
    const NSRDSfunc2& vapourThermalConductivity,
    const NSRDSfunc0& surfaceTension,
    const APIdiffCoefFunc& vapourDiffussivity
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
    K_(thermalConductivity),
    Kg_(vapourThermalConductivity),
    sigma_(surfaceTension),
    D_(vapourDiffussivity)
{}


Foam::C2H5OH::C2H5OH(Istream& is)
:
    liquidProperties(is),
    rho_(is),
    pv_(is),
    hl_(is),
    Cp_(is),
    h_(is),
    Cpg_(is),
    B_(is),
    mu_(is),
    mug_(is),
    K_(is),
    Kg_(is),
    sigma_(is),
    D_(is)
{}


Foam::C2H5OH::C2H5OH(const dictionary& dict)
:
    liquidProperties(dict),
    rho_(dict.subDict("rho")),
    pv_(dict.subDict("pv")),
    hl_(dict.subDict("hl")),
    Cp_(dict.subDict("Cp")),
    h_(dict.subDict("h")),
    Cpg_(dict.subDict("Cpg")),
    B_(dict.subDict("B")),
    mu_(dict.subDict("mu")),
    mug_(dict.subDict("mug")),
    K_(dict.subDict("K")),
    Kg_(dict.subDict("Kg")),
    sigma_(dict.subDict("sigma")),
    D_(dict.subDict("D"))
{}


Foam::C2H5OH::C2H5OH(const C2H5OH& liq)
:
    liquidProperties(liq),
    rho_(liq.rho_),
    pv_(liq.pv_),
    hl_(liq.hl_),
    Cp_(liq.Cp_),
    h_(liq.h_),
    Cpg_(liq.Cpg_),
    B_(liq.B_),
    mu_(liq.mu_),
    mug_(liq.mug_),
    K_(liq.K_),
    Kg_(liq.Kg_),
    sigma_(liq.sigma_),
    D_(liq.D_)
{}


// ************************************************************************* //
