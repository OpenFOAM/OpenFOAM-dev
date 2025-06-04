/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024-2025 OpenFOAM Foundation
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

#include "ParkRogakGroupSurfaceAreaVolumeRatioFvScalarFieldSource.H"
#include "fractal.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::ParkRogakGroupSurfaceAreaVolumeRatioFvScalarFieldSource::value
(
    const label j,
    const fvSource& source
) const
{
    const populationBalanceModel& popBal = this->popBal();
    const label i = this->i();

    const populationBalance::shapeModels::fractal& fractal =
        refCast<const populationBalance::shapeModels::fractal>(popBal.shape());

    const volScalarField::Internal& kappaj = fractal.fld(j);

    const dimensionedScalar& vi = popBal.v(i);
    const dimensionedScalar& vj = popBal.v(j);

    const volScalarField::Internal a(kappaj*vj);

    const volScalarField::Internal da1
    (
        (2.0/3.0)
       *(vi - vj)
       *(kappaj + fractal.Df(j)*(1/fractal.d(j)()() - kappaj/6))
    );

    const volScalarField::Internal dp
    (
        6/kappaj + 6*((vi - vj)*a - vj*da1)/sqr(a)
    );

    tmp<volScalarField::Internal> np
    (
        6*vi/constant::mathematical::pi/pow3(dp)
    );

    tmp<volScalarField::Internal> dc
    (
        dp*pow(np/fractal.alphaC(i), 1/fractal.Df(i))
    );

    tmp<volScalarField::Internal> da2
    (
        (vi - vj)*(4/dp + 2*fractal.Df(i)/3*(1/dc - 1/dp))
    );

    return (a + 0.5*da1 + 0.5*da2)/vi;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeTypeFieldSource
    (
        fvScalarFieldSource,
        ParkRogakGroupSurfaceAreaVolumeRatioFvScalarFieldSource
    );

    // Backwards compatible lookup as "ParkRogakSurfaceAreaVolumeRatio"
    addNamedToRunTimeSelectionTable
    (
        fvScalarFieldSource,
        ParkRogakGroupSurfaceAreaVolumeRatioFvScalarFieldSource,
        dictionary,
        ParkRogakSurfaceAreaVolumeRatio
    );
}

// ************************************************************************* //
