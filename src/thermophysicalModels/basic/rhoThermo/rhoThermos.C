/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "rhoThermo.H"
#include "makeThermo.H"

#include "specie.H"
#include "perfectGas.H"
#include "incompressiblePerfectGas.H"
#include "rhoConst.H"
#include "perfectFluid.H"
#include "PengRobinsonGas.H"
#include "adiabaticPerfectFluid.H"

#include "hConstThermo.H"
#include "janafThermo.H"
#include "sensibleEnthalpy.H"
#include "sensibleInternalEnergy.H"
#include "thermo.H"

#include "constTransport.H"
#include "sutherlandTransport.H"

#include "icoPolynomial.H"
#include "hPolynomialThermo.H"
#include "polynomialTransport.H"

#include "heRhoThermo.H"
#include "pureMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    janafThermo,
    perfectGas,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleEnthalpy,
    hConstThermo,
    rhoConst,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectFluid,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleEnthalpy,
    hConstThermo,
    adiabaticPerfectFluid,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    polynomialTransport,
    sensibleEnthalpy,
    hPolynomialThermo,
    icoPolynomial,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleEnthalpy,
    hConstThermo,
    incompressiblePerfectGas,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    hConstThermo,
    incompressiblePerfectGas,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    janafThermo,
    incompressiblePerfectGas,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    hConstThermo,
    PengRobinsonGas,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    polynomialTransport,
    sensibleEnthalpy,
    hPolynomialThermo,
    PengRobinsonGas,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    polynomialTransport,
    sensibleEnthalpy,
    janafThermo,
    PengRobinsonGas,
    specie
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleInternalEnergy,
    hConstThermo,
    perfectGas,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    sutherlandTransport,
    sensibleInternalEnergy,
    hConstThermo,
    perfectGas,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    sutherlandTransport,
    sensibleInternalEnergy,
    janafThermo,
    perfectGas,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleInternalEnergy,
    hConstThermo,
    rhoConst,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleInternalEnergy,
    hConstThermo,
    perfectFluid,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleInternalEnergy,
    hConstThermo,
    adiabaticPerfectFluid,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    polynomialTransport,
    sensibleInternalEnergy,
    hPolynomialThermo,
    icoPolynomial,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleInternalEnergy,
    hConstThermo,
    incompressiblePerfectGas,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    sutherlandTransport,
    sensibleInternalEnergy,
    hConstThermo,
    incompressiblePerfectGas,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    sutherlandTransport,
    sensibleInternalEnergy,
    janafThermo,
    incompressiblePerfectGas,
    specie
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
