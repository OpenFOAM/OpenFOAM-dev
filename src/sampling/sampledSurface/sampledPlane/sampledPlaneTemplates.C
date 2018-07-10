/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "sampledPlane.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledPlane::sampleField
(
    const GeometricField<Type, fvPatchField, volMesh>& vField
) const
{
    return tmp<Field<Type>>(new Field<Type>(vField, meshCells()));
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledPlane::interpolateField
(
    const interpolation<Type>& interpolator
) const
{
    // One value per point
    tmp<Field<Type>> tvalues(new Field<Type>(points().size()));
    Field<Type>& values = tvalues.ref();

    boolList pointDone(points().size(), false);

    forAll(faces(), cutFacei)
    {
        const face& f = faces()[cutFacei];

        forAll(f, faceVertI)
        {
            label pointi = f[faceVertI];

            if (!pointDone[pointi])
            {
                values[pointi] = interpolator.interpolate
                (
                    points()[pointi],
                    meshCells()[cutFacei]
                );
                pointDone[pointi] = true;
            }
        }
    }

    return tvalues;
}


// ************************************************************************* //
