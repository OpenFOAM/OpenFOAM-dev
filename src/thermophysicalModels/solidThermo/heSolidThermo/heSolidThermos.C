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

#include "makeHeSolidThermo.H"
#include "solidPressureThermo.H"
#include "heSolidThermo.H"

#include "specie.H"
#include "rhoConst.H"
#include "hConstThermo.H"
#include "hPowerThermo.H"
#include "hPolynomialThermo.H"
#include "constIsoSolidTransport.H"
#include "constAnIsoSolidTransport.H"
#include "exponentialSolidTransport.H"
#include "polynomialSolidTransport.H"
#include "pureMixture.H"
#include "sensibleEnthalpy.H"
#include "sensibleInternalEnergy.H"
#include "thermo.H"

#include "solidThermoPhysicsTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * * * Enthalpy-based * * * * * * * * * * * * * */

makeSolidThermo
(
    solidPressureThermo,
    heSolidThermo,
    pureMixture,
    constIsoSolidTransport,
    sensibleEnthalpy,
    hConstThermo,
    rhoConst,
    specie
);

makeSolidThermo
(
    solidPressureThermo,
    heSolidThermo,
    pureMixture,
    constAnIsoSolidTransport,
    sensibleEnthalpy,
    hConstThermo,
    rhoConst,
    specie
);

makeSolidThermo
(
    solidPressureThermo,
    heSolidThermo,
    pureMixture,
    exponentialSolidTransport,
    sensibleEnthalpy,
    hPowerThermo,
    rhoConst,
    specie
);

makeSolidThermoPhysicsType
(
    solidPressureThermo,
    heSolidThermo,
    pureMixture,
    hTransportThermoPoly8SolidThermoPhysics
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
