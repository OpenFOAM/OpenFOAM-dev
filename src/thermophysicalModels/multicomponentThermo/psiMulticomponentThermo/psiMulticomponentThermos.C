/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "coefficientMulticomponentMixture.H"
#include "coefficientWilkeMulticomponentMixture.H"
#include "singleComponentMixture.H"

#include "psiThermo.H"
#include "psiMulticomponentThermo.H"
#include "hePsiThermo.H"

#include "forGases.H"
#include "makeMulticomponentThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makePsiMulticomponentThermos(Mixture, ThermoPhysics)                   \
    makeMulticomponentThermos                                                  \
    (                                                                          \
        psiThermo,                                                             \
        psiMulticomponentThermo,                                               \
        hePsiThermo,                                                           \
        Mixture,                                                               \
        ThermoPhysics                                                          \
    )

#define makePsiMulticomponentThermo(Mixture, ThermoPhysics)                    \
    makeMulticomponentThermo                                                   \
    (                                                                          \
        psiMulticomponentThermo,                                               \
        hePsiThermo,                                                           \
        Mixture,                                                               \
        ThermoPhysics                                                          \
    )

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    forCoeffGases
    (
        makePsiMulticomponentThermos,
        coefficientMulticomponentMixture
    );

    forCoeffGases
    (
        makePsiMulticomponentThermos,
        coefficientWilkeMulticomponentMixture
    );

    forGases(makePsiMulticomponentThermo, singleComponentMixture);
}

// ************************************************************************* //
