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

#include "sampledThresholdCellFaces.H"

#include "thresholdCellFaces.H"
#include "volFieldsFwd.H"
#include "pointFields.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::sampledThresholdCellFaces::sampleField
(
    const GeometricField<Type, fvPatchField, volMesh>& vField
) const
{
    // Recreate geometry if time has changed
    updateGeometry();

    return tmp<Field<Type> >(new Field<Type>(vField, meshCells_));
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::sampledThresholdCellFaces::interpolateField
(
    const interpolation<Type>& interpolator
) const
{
    // Recreate geometry if time has changed
    updateGeometry();

    // One value per point
    tmp<Field<Type> > tvalues(new Field<Type>(points().size()));
    Field<Type>& values = tvalues();

    boolList pointDone(points().size(), false);

    forAll(faces(), cutFaceI)
    {
        const face& f = faces()[cutFaceI];

        forAll(f, faceVertI)
        {
            label pointI = f[faceVertI];

            if (!pointDone[pointI])
            {
                values[pointI] = interpolator.interpolate
                (
                    points()[pointI],
                    meshCells_[cutFaceI]
                );
                pointDone[pointI] = true;
            }
        }
    }

    return tvalues;
}


// ************************************************************************* //
