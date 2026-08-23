/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "fluidMulticomponentCarrierThermo.H"

#include "coefficientMulticomponentMixture.H"
#include "coefficientWilkeMulticomponentMixture.H"
#include "valueMulticomponentMixture.H"

#include "forGases.H"
#include "forLiquids.H"
#include "forTabulated.H"

#include "makeCarrierThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    forCoeffGases
    (
        makeCarrierThermos,
        fluidCarrierThermo,
        fluidMulticomponentCarrierThermo,
        coefficientMulticomponentMixture
    );
    forCoeffGases
    (
        makeCarrierThermos,
        fluidCarrierThermo,
        fluidMulticomponentCarrierThermo,
        coefficientWilkeMulticomponentMixture
    );

    forCoeffLiquids
    (
        makeCarrierThermos,
        fluidCarrierThermo,
        fluidMulticomponentCarrierThermo,
        coefficientMulticomponentMixture
    );
    forLiquids
    (
        makeCarrierThermos,
        fluidCarrierThermo,
        fluidMulticomponentCarrierThermo,
        valueMulticomponentMixture
    );

    forTabulated
    (
        makeCarrierThermos,
        fluidCarrierThermo,
        fluidMulticomponentCarrierThermo,
        valueMulticomponentMixture
    );
}

// ************************************************************************* //
