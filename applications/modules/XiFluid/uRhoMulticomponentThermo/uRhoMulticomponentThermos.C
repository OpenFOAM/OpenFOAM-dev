/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025-2026 OpenFOAM Foundation
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

#include "uRhoMulticomponentThermo.H"

#include "uHomogeneousMixture.H"
#include "uInhomogeneousMixture.H"
#include "uInhomogeneousEGRMixture.H"
#include "uMulticomponentMixture.H"

#include "forGases.H"

#include "makeThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeURhoMulticomponentThermos(Mixture, ThermoPhysics)                  \
                                                                               \
    defineThermo(uRhoMulticomponentThermo, Mixture, ThermoPhysics);            \
                                                                               \
    addThermo(basicThermo, uRhoMulticomponentThermo, Mixture, ThermoPhysics);  \
    addThermo(fluidThermo, uRhoMulticomponentThermo, Mixture, ThermoPhysics);  \
    addThermo                                                                  \
    (                                                                          \
        rhoFluidThermo,                                                        \
        uRhoMulticomponentThermo,                                              \
        Mixture,                                                               \
        ThermoPhysics                                                          \
    );                                                                         \
    addThermo                                                                  \
    (                                                                          \
        uRhoMulticomponentThermo,                                              \
        uRhoMulticomponentThermo,                                              \
        Mixture,                                                               \
        ThermoPhysics                                                          \
    )

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    forCoeffEnthalpyGases
    (
        makeURhoMulticomponentThermos,
        uHomogeneousMixture
    );
    forCoeffEnthalpyGases
    (
        makeURhoMulticomponentThermos,
        uInhomogeneousMixture
    );
    forCoeffEnthalpyGases
    (
        makeURhoMulticomponentThermos,
        uInhomogeneousEGRMixture
    );
    forCoeffEnthalpyGases
    (
        makeURhoMulticomponentThermos,
        uMulticomponentMixture
    );
}

// ************************************************************************* //
