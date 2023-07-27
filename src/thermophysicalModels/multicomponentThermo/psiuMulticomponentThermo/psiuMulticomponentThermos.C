/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "egrMixture.H"
#include "homogeneousMixture.H"
#include "inhomogeneousMixture.H"
#include "veryInhomogeneousMixture.H"

#include "forAbsoluteGases.H"

#include "makeThermo.H"

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

namespace Foam
{
    forAbsoluteGases(makePsiuMulticomponentThermos, egrMixture);
    forAbsoluteGases(makePsiuMulticomponentThermos, homogeneousMixture);
    forAbsoluteGases(makePsiuMulticomponentThermos, inhomogeneousMixture);
    forAbsoluteGases(makePsiuMulticomponentThermos, veryInhomogeneousMixture);
}

// ************************************************************************* //
