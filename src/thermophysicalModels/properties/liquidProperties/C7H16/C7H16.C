/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "C7H16.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(C7H16, 0);
    addToRunTimeSelectionTable(liquidProperties, C7H16,);
    addToRunTimeSelectionTable(liquidProperties, C7H16, Istream);
    addToRunTimeSelectionTable(liquidProperties, C7H16, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::C7H16::C7H16()
:
    liquidProperties
    (
        100.204,
        540.20,
        2.74e+6,
        0.428,
        0.261,
        182.57,
        1.8269e-1,
        371.58,
        0.0,
        0.3495,
        1.52e+4
    ),
    rho_(61.38396836, 0.26211, 540.2, 0.28141),
    pv_(87.829, -6996.4, -9.8802, 7.2099e-06, 2.0),
    hl_(540.20, 499121.791545248, 0.38795, 0.0, 0.0, 0.0),
    Cp_
    (
        540.20,
        6.11976102401216,
        3137.69909384855,
        182.274175063868,
       -254.530511150515
    ),
    h_
    (
       -3.1469964e+6,
        7.3072e+3,
       -3.52884e+1,
        1.10637e-1,
       -1.634831e-4,
        9.64941e-8
    ),
    Cpg_(1199.05392998284, 3992.85457666361, 1676.6, 2734.42177956968, 756.4),
    B_
    (
        0.00274040956448844,
       -2.90407568560137,
       -440900.562851782,
       -8.78208454752305e+17,
        1.28238393676899e+20
    ),
    mu_(-24.451, 1533.1, 2.0087, 0.0, 0.0),
    mug_(6.672e-08, 0.82837, 85.752, 0.0),
    K_(0.215, -0.000303, 0.0, 0.0, 0.0, 0.0),
    Kg_(-0.070028, 0.38068, -7049.9, -2400500.0),
    sigma_(540.20, 0.054143, 1.2512, 0.0, 0.0, 0.0),
    D_(147.18, 20.1, 100.204, 28.0)
{}


Foam::C7H16::C7H16
(
    const liquidProperties& l,
    const NSRDSfunc5& density,
    const NSRDSfunc1& vapourPressure,
    const NSRDSfunc6& heatOfVapourisation,
    const NSRDSfunc14& heatCapacity,
    const NSRDSfunc0& enthalpy,
    const NSRDSfunc7& idealGasHeatCapacity,
    const NSRDSfunc4& secondVirialCoeff,
    const NSRDSfunc1& dynamicViscosity,
    const NSRDSfunc2& vapourDynamicViscosity,
    const NSRDSfunc0& thermalConductivity,
    const NSRDSfunc2& vapourThermalConductivity,
    const NSRDSfunc6& surfaceTension,
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


Foam::C7H16::C7H16(Istream& is)
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


Foam::C7H16::C7H16(const dictionary& dict)
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


Foam::C7H16::C7H16(const C7H16& liq)
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
