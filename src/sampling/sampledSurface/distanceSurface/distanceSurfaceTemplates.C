/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "distanceSurface.H"
#include "volFieldsFwd.H"
#include "pointFields.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::distanceSurface::sampleField
(
    const GeometricField<Type, fvPatchField, volMesh>& vField
) const
{
    if (cell_)
    {
        return tmp<Field<Type> >
        (
            new Field<Type>(vField, isoSurfCellPtr_().meshCells())
        );
    }
    else
    {
        return tmp<Field<Type> >
        (
            new Field<Type>(vField, isoSurfPtr_().meshCells())
        );
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::distanceSurface::interpolateField
(
    const interpolation<Type>& interpolator
) const
{
    const fvMesh& fvm = static_cast<const fvMesh&>(mesh());

    // Get fields to sample. Assume volPointInterpolation!
    const GeometricField<Type, fvPatchField, volMesh>& volFld =
        interpolator.psi();

    tmp<GeometricField<Type, pointPatchField, pointMesh> > pointFld
    (
        volPointInterpolation::New(fvm).interpolate(volFld)
    );

    // Sample.
    if (cell_)
    {
        return isoSurfCellPtr_().interpolate
        (
            (
                average_
              ? pointAverage(pointFld())()
              : volFld
            ),
            pointFld()
        );
    }
    else
    {
        return isoSurfPtr_().interpolate
        (
            (
                average_
              ? pointAverage(pointFld())()
              : volFld
            ),
            pointFld()
        );
    }
}


// ************************************************************************* //
