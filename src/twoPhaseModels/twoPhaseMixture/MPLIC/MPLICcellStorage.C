/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2021 OpenFOAM Foundation
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

#include "MPLICcellStorage.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::boolList Foam::MPLICcellStorage::calcIsOwner
(
    const primitiveMesh& mesh,
    const label celli
) const
{
    const cell& c = mesh.cells()[celli];
    boolList cOwns(c.size(), false);
    forAll(c, i)
    {
        const label ow = mesh.faceOwner()[c[i]];
        if (ow == celli)
        {
            cOwns[i] = true;
        }
    }
    return cOwns;
}


Foam::scalar Foam::MPLICcellStorage::calcAlphaMin() const
{
    // Initialise with the first value in the list
    scalar cellAlphaMin(pointsAlpha_.first());

    // Loop through the rest and compare element-wise
    for (label i = 1; i < cPoints_.size(); ++i)
    {
        const label pI = cPoints_[i];
        cellAlphaMin = min(cellAlphaMin, pointsAlpha_[pI]);
    }

    // Return minimum value
    return cellAlphaMin;
}


Foam::scalar Foam::MPLICcellStorage::calcAlphaMax() const
{
    // Initialise with the first value in the list
    scalar cellAlphaMax(pointsAlpha_.first());

    // Loop through the rest and compare element-wise
    for (label i = 1; i < cPoints_.size(); ++i)
    {
        const label pI = cPoints_[i];
        cellAlphaMax = max(cellAlphaMax, pointsAlpha_[pI]);
    }

    // Return maximum value
    return cellAlphaMax;
}


Foam::scalarField Foam::MPLICcellStorage::calcFacesAlphaMin() const
{
    scalarField facesAlphaMin(cFaces_.size());

    forAll(cFaces_, facei)
    {
        // Face
        const face& f = faces_[cFaces_[facei]];

        // Initialise with the first value in the list
        scalar fAlphaMin(pointsAlpha_[f.first()]);

        // Loop through the rest and compare element-wise
        for (label i = 1; i < f.size(); ++i)
        {
            fAlphaMin = min(fAlphaMin, pointsAlpha_[f[i]]);
        }

        facesAlphaMin[facei] = fAlphaMin;
    }

    // Return minimum value
    return facesAlphaMin;
}


Foam::scalarField Foam::MPLICcellStorage::calcFacesAlphaMax() const
{
    scalarField facesAlphaMax(cFaces_.size());

    forAll(cFaces_, facei)
    {
        // Face
        const face& f = faces_[cFaces_[facei]];

        // Initialise with the first value in the list
        scalar fAlphaMax(pointsAlpha_[f.first()]);

        // Loop through the rest and compare element-wise
        for (label i = 1; i < f.size(); ++i)
        {
            fAlphaMax = max(fAlphaMax, pointsAlpha_[f[i]]);
        }

        facesAlphaMax[facei] = fAlphaMax;
    }

    // Return maximum
    return facesAlphaMax;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MPLICcellStorage::MPLICcellStorage
(
    const primitiveMesh& mesh,
    const scalarField& pointsAlpha,
    const vectorField& pointsU,
    const scalar cellAlpha,
    const vector& cellU,
    const label celli
)
:
    points_(mesh.points()),
    faces_(mesh.faces()),
    edges_(mesh.edges()),
    edgeFaces_(mesh.faceEdges()),
    cPoints_(mesh.cellPoints()[celli]),
    cFaces_(mesh.cells()[celli]),
    cEdges_(mesh.cellEdges()[celli]),
    pointsAlpha_(pointsAlpha),
    pointsU_(pointsU),
    cellAlpha_(cellAlpha),
    celllU_(cellU),
    owns_(calcIsOwner(mesh, celli)),
    volume_(mesh.cellVolumes()[celli]),
    centre_(mesh.cellCentres()[celli]),
    Sf_(mesh.faceAreas(), cFaces_),
    Cf_(mesh.faceCentres(), cFaces_),
    magSf_(mesh.magFaceAreas(), cFaces_),
    cellAlphaMin_(calcAlphaMin()),
    cellAlphaMax_(calcAlphaMax()),
    facesAlphaMin_(calcFacesAlphaMin()),
    facesAlphaMax_(calcFacesAlphaMax()),
    faceEdges_(mesh.faceEdges())
{}


// ************************************************************************* //
