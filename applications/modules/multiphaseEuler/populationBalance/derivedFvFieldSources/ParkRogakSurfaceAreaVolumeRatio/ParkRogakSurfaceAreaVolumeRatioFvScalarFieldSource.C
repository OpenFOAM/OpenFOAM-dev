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

#include "ParkRogakSurfaceAreaVolumeRatioFvScalarFieldSource.H"
#include "fractal.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::ParkRogakSurfaceAreaVolumeRatioFvScalarFieldSource::value
(
    const label deltai,
    const fvSource& source
) const
{
    const diameterModels::sizeGroup& fi = this->fi();
    const diameterModels::sizeGroup& fj = this->fi(deltai);

    const volScalarField::Internal& kappaj =
        fld<diameterModels::shapeModels::fractal>(deltai);

    const diameterModels::shapeModels::fractal& fractali =
        model<diameterModels::shapeModels::fractal>();
    const diameterModels::shapeModels::fractal& fractalj =
        model<diameterModels::shapeModels::fractal>(deltai);

    const volScalarField::Internal a(kappaj*fj.x());

    const dimensionedScalar dv(fi.x() - fj.x());

    const volScalarField::Internal da1
    (
        (2.0/3.0)
       *dv
       *(kappaj + fractalj.Df()*(1/fractalj.d()()() - kappaj/6))
    );

    const volScalarField::Internal dp
    (
        6/kappaj + 6*(dv*a - fj.x()*da1)/sqr(a)
    );

    tmp<volScalarField::Internal> np
    (
        6*fi.x()/constant::mathematical::pi/pow3(dp)
    );

    tmp<volScalarField::Internal> dc
    (
        dp*pow(np/fractali.alphaC(), 1/fractali.Df())
    );

    tmp<volScalarField::Internal> da2
    (
        dv*(4/dp + 2*fractali.Df()/3*(1/dc - 1/dp))
    );

    return (a + 0.5*da1 + 0.5*da2)/fi.x();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeTypeFieldSource
    (
        fvScalarFieldSource,
        ParkRogakSurfaceAreaVolumeRatioFvScalarFieldSource
    );
}

// ************************************************************************* //
