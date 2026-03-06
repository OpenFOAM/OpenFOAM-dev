/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024-2026 OpenFOAM Foundation
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

#include "hardSphereGroupSurfaceAreaVolumeRatioFvScalarFieldSource.H"
#include "fractal.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::fvMesh>>
Foam::hardSphereGroupSurfaceAreaVolumeRatioFvScalarFieldSource::value
(
    const label j,
    const fvSource& source
) const
{
    const populationBalanceModel& popBal = this->popBal();

    const populationBalance::shapeModels::fractal& fractal =
        refCast<const populationBalance::shapeModels::fractal>(popBal.shape());

    const volScalarField::Internal& kappaj = fractal.fld(j);

    return kappaj*popBal.dSph(j)/popBal.dSph(i());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeTypeFieldSource
    (
        fvScalarFieldSource,
        hardSphereGroupSurfaceAreaVolumeRatioFvScalarFieldSource
    );

    // Backwards compatible lookup as "hardSphereSurfaceAreaVolumeRatio"
    addNamedToRunTimeSelectionTable
    (
        fvScalarFieldSource,
        hardSphereGroupSurfaceAreaVolumeRatioFvScalarFieldSource,
        dictionary,
        hardSphereSurfaceAreaVolumeRatio
    );
}

// ************************************************************************* //
