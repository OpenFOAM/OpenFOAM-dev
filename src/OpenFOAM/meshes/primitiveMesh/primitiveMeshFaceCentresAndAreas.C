/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

Description
    Calculate the face centres and areas.

    Calculate the centre by breaking the face into triangles using the face
    centre and area-weighted averaging their centres.  This method copes with
    small face-concavity.

\*---------------------------------------------------------------------------*/

#include "primitiveMesh.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::primitiveMesh::calcFaceCentresAndAreas() const
{
    if (debug)
    {
        Pout<< "primitiveMesh::calcFaceCentresAndAreas() : "
            << "Calculating face centres and face areas"
            << endl;
    }

    // It is an error to attempt to recalculate faceCentres
    // if the pointer is already set
    if (faceCentresPtr_ || faceAreasPtr_ || magFaceAreasPtr_)
    {
        FatalErrorInFunction
            << "Face centres or face areas already calculated"
            << abort(FatalError);
    }

    faceCentresPtr_ = new vectorField(nFaces());
    vectorField& fCtrs = *faceCentresPtr_;

    faceAreasPtr_ = new vectorField(nFaces());
    vectorField& fAreas = *faceAreasPtr_;

    magFaceAreasPtr_ = new scalarField(nFaces());
    scalarField& magfAreas = *magFaceAreasPtr_;

    makeFaceCentresAndAreas(points(), fCtrs, fAreas, magfAreas);

    if (debug)
    {
        Pout<< "primitiveMesh::calcFaceCentresAndAreas() : "
            << "Finished calculating face centres and face areas"
            << endl;
    }
}


void Foam::primitiveMesh::makeFaceCentresAndAreas
(
    const pointField& p,
    vectorField& fCtrs,
    vectorField& fAreas,
    scalarField& magfAreas
) const
{
    const faceList& fs = faces();

    forAll(fs, facei)
    {
        const Tuple2<vector, point> areaAndCentre =
            face::areaAndCentre(UIndirectList<point>(p, fs[facei]));

        fCtrs[facei] = areaAndCentre.second();
        fAreas[facei] = areaAndCentre.first();
        magfAreas[facei] = max(mag(fAreas[facei]), rootVSmall);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::vectorField& Foam::primitiveMesh::faceCentres() const
{
    if (!faceCentresPtr_)
    {
        calcFaceCentresAndAreas();
    }

    return *faceCentresPtr_;
}


const Foam::vectorField& Foam::primitiveMesh::faceAreas() const
{
    if (!faceAreasPtr_)
    {
        calcFaceCentresAndAreas();
    }

    return *faceAreasPtr_;
}


const Foam::scalarField& Foam::primitiveMesh::magFaceAreas() const
{
    if (!magFaceAreasPtr_)
    {
        calcFaceCentresAndAreas();
    }

    return *magFaceAreasPtr_;
}


// ************************************************************************* //
