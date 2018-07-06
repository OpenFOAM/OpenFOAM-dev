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

#include "sampledTriSurfaceMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledTriSurfaceMesh::sampleField
(
    const GeometricField<Type, fvPatchField, volMesh>& vField
) const
{
    // One value per face
    tmp<Field<Type>> tvalues(new Field<Type>(sampleElements_.size()));
    Field<Type>& values = tvalues.ref();

    if (sampleSource_ == cells || sampleSource_ == insideCells)
    {
        // Sample cells

        forAll(sampleElements_, triI)
        {
            values[triI] = vField[sampleElements_[triI]];
        }
    }
    else
    {
        // Sample boundary faces

        const polyBoundaryMesh& pbm = mesh().boundaryMesh();
        label nBnd = mesh().nFaces()-mesh().nInternalFaces();

        // Create flat boundary field

        Field<Type> bVals(nBnd, Zero);

        forAll(vField.boundaryField(), patchi)
        {
            label bFacei = pbm[patchi].start() - mesh().nInternalFaces();

            SubList<Type>
            (
                bVals,
                vField.boundaryField()[patchi].size(),
                bFacei
            ) = vField.boundaryField()[patchi];
        }

        // Sample in flat boundary field

        forAll(sampleElements_, triI)
        {
            label facei = sampleElements_[triI];
            values[triI] = bVals[facei-mesh().nInternalFaces()];
        }
    }

    return tvalues;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::sampledTriSurfaceMesh::interpolateField
(
    const interpolation<Type>& interpolator
) const
{
    // One value per vertex
    tmp<Field<Type>> tvalues(new Field<Type>(sampleElements_.size()));
    Field<Type>& values = tvalues.ref();

    if (sampleSource_ == cells || sampleSource_ == insideCells)
    {
        // Sample cells.

        forAll(sampleElements_, pointi)
        {
            values[pointi] = interpolator.interpolate
            (
                samplePoints_[pointi],
                sampleElements_[pointi]
            );
        }
    }
    else
    {
        // Sample boundary faces.

        forAll(samplePoints_, pointi)
        {
            label facei = sampleElements_[pointi];

            values[pointi] = interpolator.interpolate
            (
                samplePoints_[pointi],
                mesh().faceOwner()[facei],
                facei
            );
        }
    }

    return tvalues;
}


// ************************************************************************* //
