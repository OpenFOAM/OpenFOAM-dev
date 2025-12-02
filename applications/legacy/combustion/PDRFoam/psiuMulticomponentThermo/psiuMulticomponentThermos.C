/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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

#include "psiuMulticomponentThermo.H"

#include "homogeneousMixture.H"
#include "leanInhomogeneousMixture.H"
#include "inhomogeneousMixture.H"
#include "inhomogeneousEGRMixture.H"

#include "specie.H"

#include "perfectGas.H"

#include "hConstThermo.H"
#include "janafThermo.H"

#include "absoluteEnthalpy.H"

#include "constTransport.H"
#include "sutherlandTransport.H"

#include "thermo.H"

#include "forThermo.H"

#include "makeThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define forAbsoluteEnthalpyGasEqns(Mu, He, Cp, Macro, Args...)                 \
    forThermo(Mu, He, Cp, perfectGas, specie, Macro, Args)

#define forAbsoluteEnthalpyGasEnergiesAndThermos(Mu, Args...)                  \
    forAbsoluteEnthalpyGasEqns(Mu, absoluteEnthalpy, hConstThermo, Args);      \
    forAbsoluteEnthalpyGasEqns(Mu, absoluteEnthalpy, janafThermo, Args)

#define forAbsoluteEnthalpyGasTransports(Macro, Args...)                       \
    forAbsoluteEnthalpyGasEnergiesAndThermos(constTransport, Macro, Args);     \
    forAbsoluteEnthalpyGasEnergiesAndThermos(sutherlandTransport, Macro, Args)

#define forAbsoluteEnthalpyGases(Macro, Args...)                               \
    forAbsoluteEnthalpyGasTransports(Macro, Args)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makePsiuMulticomponentThermos(Mixture, ThermoPhysics)                  \
                                                                               \
    defineThermo(psiuMulticomponentThermo, Mixture, ThermoPhysics);            \
                                                                               \
    addThermo(basicThermo, psiuMulticomponentThermo, Mixture, ThermoPhysics);  \
    addThermo(fluidThermo, psiuMulticomponentThermo, Mixture, ThermoPhysics);  \
    addThermo(psiThermo, psiuMulticomponentThermo, Mixture, ThermoPhysics);    \
    addThermo                                                                  \
    (                                                                          \
        psiuMulticomponentThermo,                                              \
        psiuMulticomponentThermo,                                              \
        Mixture,                                                               \
        ThermoPhysics                                                          \
    )

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    forAbsoluteEnthalpyGases
    (
        makePsiuMulticomponentThermos,
        homogeneousMixture
    );
    forAbsoluteEnthalpyGases
    (
        makePsiuMulticomponentThermos,
        leanInhomogeneousMixture
    );
    forAbsoluteEnthalpyGases
    (
        makePsiuMulticomponentThermos,
        inhomogeneousMixture
    );
    forAbsoluteEnthalpyGases
    (
        makePsiuMulticomponentThermos,
        inhomogeneousEGRMixture
    );
}

// ************************************************************************* //
