/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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

#include "sampledIsoSurfaceSurface.H"
#include "interpolationVolPointInterpolation.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledSurfaces::sampledIsoSurfaceSurface::sampleField
(
    const VolField<Type>& vField
) const
{
    update();

    return isoSurfPtr_->sample(vField.primitiveField());
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledSurfaces::sampledIsoSurfaceSurface::interpolateField
(
    const interpolation<Type>& interpolator
) const
{
    update();

    if (isA<interpolationVolPointInterpolation<Type>>(interpolator))
    {
        const PointField<Type>& pField =
            refCast<const interpolationVolPointInterpolation<Type>>
            (interpolator).psip();

        return isoSurfPtr_().interpolate(pField.primitiveField());
    }
    else
    {
        const pointField& points = isoSurfPtr_->points();
        const faceList& faces = isoSurfPtr_->faces();
        const labelList& faceCells = isoSurfPtr_->faceCells();

        labelList pointCells(points.size());
        forAll(faces, facei)
        {
            forAll(faces[facei], facePointi)
            {
                pointCells[faces[facei][facePointi]] = faceCells[facei];
            }
        }

        return interpolator.interpolate(points, pointCells);
    }
}


// ************************************************************************* //
