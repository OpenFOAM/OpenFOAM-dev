/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "ubRhoMulticomponentThermo.H"

#include "uHomogeneousMixture.H"
#include "bHomogeneousMixture.H"
#include "uInhomogeneousMixture.H"
#include "bInhomogeneousMixture.H"
#include "uInhomogeneousEGRMixture.H"

#include "forGases.H"

#include "makeThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeUBRhoMulticomponentThermos(Mixture, ThermoPhysics)                 \
                                                                               \
    defineThermo(ubRhoMulticomponentThermo, Mixture, ThermoPhysics);           \
                                                                               \
    addThermo(basicThermo, ubRhoMulticomponentThermo, Mixture, ThermoPhysics); \
    addThermo(fluidThermo, ubRhoMulticomponentThermo, Mixture, ThermoPhysics); \
    addThermo                                                                  \
    (                                                                          \
        rhoFluidThermo,                                                        \
        ubRhoMulticomponentThermo,                                             \
        Mixture,                                                               \
        ThermoPhysics                                                          \
    );                                                                         \
    addThermo                                                                  \
    (                                                                          \
        ubRhoMulticomponentThermo,                                             \
        ubRhoMulticomponentThermo,                                             \
        Mixture,                                                               \
        ThermoPhysics                                                          \
    )

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    forCoeffEnthalpyGases
    (
        makeUBRhoMulticomponentThermos,
        uHomogeneousMixture
    );
    forCoeffEnthalpyGases
    (
        makeUBRhoMulticomponentThermos,
        bHomogeneousMixture
    );
    forCoeffEnthalpyGases
    (
        makeUBRhoMulticomponentThermos,
        uInhomogeneousMixture
    );
    forCoeffEnthalpyGases
    (
        makeUBRhoMulticomponentThermos,
        bInhomogeneousMixture
    );
    forCoeffEnthalpyGases
    (
        makeUBRhoMulticomponentThermos,
        uInhomogeneousEGRMixture
    );
}

// ************************************************************************* //
