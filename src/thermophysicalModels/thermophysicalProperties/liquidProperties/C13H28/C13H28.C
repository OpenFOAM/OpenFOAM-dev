/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "C13H28.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(C13H28, 0);
    addToRunTimeSelectionTable(liquidProperties, C13H28,);
    addToRunTimeSelectionTable(liquidProperties, C13H28, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::C13H28::C13H28()
:
    liquidProperties
    (
        184.365,
        675.80,
        1.7225e+6,
        0.77,
        0.236,
        267.76,
        3.801e-1,
        508.62,
        0.0,
        0.6186,
        1.5901e+4
    ),
    rho_(59.513022, 0.2504, 675.8, 0.312),
    pv_(118.27, -11432, -13.769, 5.9641e-06, 2.0),
    hl_(675.80, 444227.48352453, 0.4162, 0.0, 0.0, 0.0),
    Cp_
    (
        4275.05220622135,
       -16.6539202126217,
        0.0325755973205326,
        0.0,
        0.0,
        0.0
    ),
    h_
    (
       -2860442.0545124,
        4275.05220622135,
       -8.32696010631085,
        0.0108585324401775,
        0.0,
        0.0
    ),
    Cpg_(1136.87522035093, 3641.14663846175, -1443, 2277.00485450058, -683.0),
    B_
    (
        0.00246321156401703,
       -2.66601578390692,
       -1249532.17801643,
       -1.0460770753668e+19,
        1.90117430097904e+21
    ),
    mu_(-23.341, 2121.9, 1.7208, 0.0, 0.0),
    mug_(3.5585e-08, 0.8987, 165.3, 0.0),
    kappa_(0.1981, -0.0002046, 0.0, 0.0, 0.0, 0.0),
    kappag_(5.3701e-06, 1.4751, 599.09, 0.0),
    sigma_(675.80, 0.05561, 1.3361, 0.0, 0.0, 0.0),
    D_(147.18, 20.1, 184.365, 28.0) // note: Same as nHeptane
{}


Foam::C13H28::C13H28
(
    const liquidProperties& l,
    const thermophysicalFunctions::NSRDS5& density,
    const thermophysicalFunctions::NSRDS1& vapourPressure,
    const thermophysicalFunctions::NSRDS6& heatOfVapourisation,
    const thermophysicalFunctions::NSRDS0& heatCapacity,
    const thermophysicalFunctions::NSRDS0& enthalpy,
    const thermophysicalFunctions::NSRDS7& idealGasHeatCapacity,
    const thermophysicalFunctions::NSRDS4& secondVirialCoeff,
    const thermophysicalFunctions::NSRDS1& dynamicViscosity,
    const thermophysicalFunctions::NSRDS2& vapourDynamicViscosity,
    const thermophysicalFunctions::NSRDS0& thermalConductivity,
    const thermophysicalFunctions::NSRDS2& vapourThermalConductivity,
    const thermophysicalFunctions::NSRDS6& surfaceTension,
    const thermophysicalFunctions::APIdiffCoef& vapourDiffusivity
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


Foam::C13H28::C13H28(const dictionary& dict)
:
    C13H28()
{
    readIfPresent(*this, dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::C13H28::write(Ostream& os) const
{
    liquidProperties::write(*this, os);
}


// ************************************************************************* //
