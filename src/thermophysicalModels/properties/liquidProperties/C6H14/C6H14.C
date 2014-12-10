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

#include "C6H14.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(C6H14, 0);
    addToRunTimeSelectionTable(liquidProperties, C6H14,);
    addToRunTimeSelectionTable(liquidProperties, C6H14, Istream);
    addToRunTimeSelectionTable(liquidProperties, C6H14, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::C6H14::C6H14()
:
    liquidProperties
    (
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
    rho_(61.03399848, 0.26411, 507.6, 0.27537),
    pv_(104.65, -6995.5, -12.702, 1.2381e-05, 2.0),
    hl_(507.60, 527286.863084118, 0.39002, 0.0, 0.0, 0.0),
    Cp_
    (
        1997.28465831951,
       -2.13258758137322,
        0.0102964828202421,
        0.0,
        0.0,
        0.0
    ),
    h_
    (
       -2902186.5403246,
        1997.28465831951,
       -1.06629379068661,
        0.00343216094008069,
        0.0,
        0.0
    ),
    Cpg_(1211.4601343746, 4088.0977523005, 1694.6, 2748.99335089409, 761.6),
    B_
    (
        0.0022859927822969,
       -2.32080485512376,
       -430509.300625457,
        1.93787205402834e+17,
       -7.17128700233241e+19
    ),
    mu_(-20.715, 1207.5, 1.4993, 0.0, 0.0),
    mug_(1.7514e-07, 0.70737, 157.14, 0.0),
    K_(0.22492, -0.0003533, 0.0, 0.0, 0.0, 0.0),
    Kg_(-650.5, 0.8053, -1412100000, 0.0),
    sigma_(507.60, 0.055003, 1.2674, 0.0, 0.0, 0.0),
    D_(147.18, 20.1, 86.177, 28) // note: Same as nHeptane
{}


Foam::C6H14::C6H14
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


Foam::C6H14::C6H14(Istream& is)
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


Foam::C6H14::C6H14(const dictionary& dict)
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


Foam::C6H14::C6H14(const C6H14& liq)
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
