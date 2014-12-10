/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Calulate the face centres and areas.

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
    if (faceCentresPtr_ || faceAreasPtr_)
    {
        FatalErrorIn("primitiveMesh::calcFaceCentresAndAreas() const")
            << "Face centres or face areas already calculated"
            << abort(FatalError);
    }

    faceCentresPtr_ = new vectorField(nFaces());
    vectorField& fCtrs = *faceCentresPtr_;

    faceAreasPtr_ = new vectorField(nFaces());
    vectorField& fAreas = *faceAreasPtr_;

    makeFaceCentresAndAreas(points(), fCtrs, fAreas);

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
    vectorField& fAreas
) const
{
    const faceList& fs = faces();

    forAll(fs, facei)
    {
        const labelList& f = fs[facei];
        label nPoints = f.size();

        // If the face is a triangle, do a direct calculation for efficiency
        // and to avoid round-off error-related problems
        if (nPoints == 3)
        {
            fCtrs[facei] = (1.0/3.0)*(p[f[0]] + p[f[1]] + p[f[2]]);
            fAreas[facei] = 0.5*((p[f[1]] - p[f[0]])^(p[f[2]] - p[f[0]]));
        }
        else
        {
            vector sumN = vector::zero;
            scalar sumA = 0.0;
            vector sumAc = vector::zero;

            point fCentre = p[f[0]];
            for (label pi = 1; pi < nPoints; pi++)
            {
                fCentre += p[f[pi]];
            }

            fCentre /= nPoints;

            for (label pi = 0; pi < nPoints; pi++)
            {
                const point& nextPoint = p[f[(pi + 1) % nPoints]];

                vector c = p[f[pi]] + nextPoint + fCentre;
                vector n = (nextPoint - p[f[pi]])^(fCentre - p[f[pi]]);
                scalar a = mag(n);

                sumN += n;
                sumA += a;
                sumAc += a*c;
            }

            // This is to deal with zero-area faces. Mark very small faces
            // to be detected in e.g., processorPolyPatch.
            if (sumA < ROOTVSMALL)
            {
                fCtrs[facei] = fCentre;
                fAreas[facei] = vector::zero;
            }
            else
            {
                fCtrs[facei] = (1.0/3.0)*sumAc/sumA;
                fAreas[facei] = 0.5*sumN;
            }
        }
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


// ************************************************************************* //
