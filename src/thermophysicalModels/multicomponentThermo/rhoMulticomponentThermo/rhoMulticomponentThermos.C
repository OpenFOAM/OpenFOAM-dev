/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2023 OpenFOAM Foundation
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
#include "valueMulticomponentMixture.H"
#include "singleComponentMixture.H"

#include "rhoThermo.H"
#include "rhoMulticomponentThermo.H"
#include "heRhoThermo.H"
#include "heMulticomponentThermo.H"
#include "heFluidMulticomponentThermo.H"

#include "forGases.H"
#include "forLiquids.H"
#include "forTabulated.H"
#include "makeFluidMulticomponentThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    forCoeffGases
    (
        makeFluidMulticomponentThermos,
        rhoThermo,
        rhoMulticomponentThermo,
        coefficientMulticomponentMixture
    );
    forCoeffGases
    (
        makeFluidMulticomponentThermos,
        rhoThermo,
        rhoMulticomponentThermo,
        coefficientWilkeMulticomponentMixture
    );
    forGases
    (
        makeFluidMulticomponentThermo,
        rhoMulticomponentThermo,
        singleComponentMixture
    );

    forCoeffLiquids
    (
        makeFluidMulticomponentThermos,
        rhoThermo,
        rhoMulticomponentThermo,
        coefficientMulticomponentMixture
    );
    forLiquids
    (
        makeFluidMulticomponentThermos,
        rhoThermo,
        rhoMulticomponentThermo,
        valueMulticomponentMixture
    );
    forLiquids
    (
        makeFluidMulticomponentThermo,
        rhoMulticomponentThermo,
        singleComponentMixture
    );

    forTabulated
    (
        makeFluidMulticomponentThermos,
        rhoThermo,
        rhoMulticomponentThermo,
        valueMulticomponentMixture
    );
    forTabulated
    (
        makeFluidMulticomponentThermo,
        rhoMulticomponentThermo,
        singleComponentMixture
    );
}

// ************************************************************************* //
