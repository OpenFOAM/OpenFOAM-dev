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

#include "fluidMulticomponentSurfaceThermo.H"

#include "coefficientMulticomponentMixture.H"
#include "coefficientWilkeMulticomponentMixture.H"
#include "valueMulticomponentMixture.H"

#include "forGases.H"
#include "forLiquids.H"
#include "forTabulated.H"

#include "makeSurfaceThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    forCoeffGases
    (
        makeSurfaceThermos,
        fluidSurfaceThermo,
        fluidMulticomponentSurfaceThermo,
        coefficientMulticomponentMixture
    );
    forCoeffGases
    (
        makeSurfaceThermos,
        fluidSurfaceThermo,
        fluidMulticomponentSurfaceThermo,
        coefficientWilkeMulticomponentMixture
    );

    forCoeffLiquids
    (
        makeSurfaceThermos,
        fluidSurfaceThermo,
        fluidMulticomponentSurfaceThermo,
        coefficientMulticomponentMixture
    );
    forLiquids
    (
        makeSurfaceThermos,
        fluidSurfaceThermo,
        fluidMulticomponentSurfaceThermo,
        valueMulticomponentMixture
    );

    forTabulated
    (
        makeSurfaceThermos,
        fluidSurfaceThermo,
        fluidMulticomponentSurfaceThermo,
        valueMulticomponentMixture
    );
}

// ************************************************************************* //
