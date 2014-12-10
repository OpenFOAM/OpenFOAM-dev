/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2014 OpenFOAM Foundation
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

#include "primitiveMeshTools.H"
#include "syncTools.H"
#include "pyramidPointFaceRef.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::primitiveMeshTools::faceSkewness
(
    const primitiveMesh& mesh,
    const pointField& p,
    const vectorField& fCtrs,
    const vectorField& fAreas,

    const label faceI,
    const point& ownCc,
    const point& neiCc
)
{
    vector Cpf = fCtrs[faceI] - ownCc;
    vector d = neiCc - ownCc;

    // Skewness vector
    vector sv =
        Cpf
      - ((fAreas[faceI] & Cpf)/((fAreas[faceI] & d) + SMALL))*d;
    vector svHat = sv/(mag(sv) + VSMALL);

    // Normalisation distance calculated as the approximate distance
    // from the face centre to the edge of the face in the direction
    // of the skewness
    scalar fd = 0.2*mag(d) + VSMALL;
    const face& f = mesh.faces()[faceI];
    forAll(f, pi)
    {
        fd = max(fd, mag(svHat & (p[f[pi]] - fCtrs[faceI])));
    }

    // Normalised skewness
    return mag(sv)/fd;
}

Foam::scalar Foam::primitiveMeshTools::boundaryFaceSkewness
(
    const primitiveMesh& mesh,
    const pointField& p,
    const vectorField& fCtrs,
    const vectorField& fAreas,

    const label faceI,
    const point& ownCc
)
{
    vector Cpf = fCtrs[faceI] - ownCc;

    vector normal = fAreas[faceI];
    normal /= mag(normal) + VSMALL;
    vector d = normal*(normal & Cpf);


    // Skewness vector
    vector sv =
        Cpf
      - ((fAreas[faceI] & Cpf)/((fAreas[faceI] & d) + VSMALL))*d;
    vector svHat = sv/(mag(sv) + VSMALL);

    // Normalisation distance calculated as the approximate distance
    // from the face centre to the edge of the face in the direction
    // of the skewness
    scalar fd = 0.4*mag(d) + VSMALL;
    const face& f = mesh.faces()[faceI];
    forAll(f, pi)
    {
        fd = max(fd, mag(svHat & (p[f[pi]] - fCtrs[faceI])));
    }

    // Normalised skewness
    return mag(sv)/fd;
}


Foam::scalar Foam::primitiveMeshTools::faceOrthogonality
(
    const point& ownCc,
    const point& neiCc,
    const vector& s
)
{
    vector d = neiCc - ownCc;

    return (d & s)/(mag(d)*mag(s) + VSMALL);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::primitiveMeshTools::faceOrthogonality
(
    const primitiveMesh& mesh,
    const vectorField& areas,
    const vectorField& cc
)
{
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();

    tmp<scalarField> tortho(new scalarField(mesh.nInternalFaces()));
    scalarField& ortho = tortho();

    // Internal faces
    forAll(nei, faceI)
    {
        ortho[faceI] = faceOrthogonality
        (
            cc[own[faceI]],
            cc[nei[faceI]],
            areas[faceI]
        );
    }

    return tortho;
}


Foam::tmp<Foam::scalarField> Foam::primitiveMeshTools::faceSkewness
(
    const primitiveMesh& mesh,
    const pointField& p,
    const vectorField& fCtrs,
    const vectorField& fAreas,
    const vectorField& cellCtrs
)
{
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();

    tmp<scalarField> tskew(new scalarField(mesh.nFaces()));
    scalarField& skew = tskew();

    forAll(nei, faceI)
    {
        skew[faceI] = faceSkewness
        (
            mesh,
            p,
            fCtrs,
            fAreas,

            faceI,
            cellCtrs[own[faceI]],
            cellCtrs[nei[faceI]]
        );
    }


    // Boundary faces: consider them to have only skewness error.
    // (i.e. treat as if mirror cell on other side)

    for (label faceI = mesh.nInternalFaces(); faceI < mesh.nFaces(); faceI++)
    {
        skew[faceI] = boundaryFaceSkewness
        (
            mesh,
            p,
            fCtrs,
            fAreas,
            faceI,
            cellCtrs[own[faceI]]
        );
    }

    return tskew;
}

void Foam::primitiveMeshTools::facePyramidVolume
(
    const primitiveMesh& mesh,
    const pointField& points,
    const vectorField& ctrs,

    scalarField& ownPyrVol,
    scalarField& neiPyrVol
)
{
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    const faceList& f = mesh.faces();

    ownPyrVol.setSize(mesh.nFaces());
    neiPyrVol.setSize(mesh.nInternalFaces());

    forAll(f, faceI)
    {
        // Create the owner pyramid
        ownPyrVol[faceI] = -pyramidPointFaceRef
        (
            f[faceI],
            ctrs[own[faceI]]
        ).mag(points);

        if (mesh.isInternalFace(faceI))
        {
            // Create the neighbour pyramid - it will have positive volume
            neiPyrVol[faceI] = pyramidPointFaceRef
            (
                f[faceI],
                ctrs[nei[faceI]]
            ).mag(points);
        }
    }
}


void Foam::primitiveMeshTools::cellClosedness
(
    const primitiveMesh& mesh,
    const Vector<label>& meshD,
    const vectorField& areas,
    const scalarField& vols,

    scalarField& openness,
    scalarField& aratio
)
{
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();

    // Loop through cell faces and sum up the face area vectors for each cell.
    // This should be zero in all vector components

    vectorField sumClosed(mesh.nCells(), vector::zero);
    vectorField sumMagClosed(mesh.nCells(), vector::zero);

    forAll(own, faceI)
    {
        // Add to owner
        sumClosed[own[faceI]] += areas[faceI];
        sumMagClosed[own[faceI]] += cmptMag(areas[faceI]);
    }

    forAll(nei, faceI)
    {
        // Subtract from neighbour
        sumClosed[nei[faceI]] -= areas[faceI];
        sumMagClosed[nei[faceI]] += cmptMag(areas[faceI]);
    }


    label nDims = 0;
    for (direction dir = 0; dir < vector::nComponents; dir++)
    {
        if (meshD[dir] == 1)
        {
            nDims++;
        }
    }


    // Check the sums
    openness.setSize(mesh.nCells());
    aratio.setSize(mesh.nCells());

    forAll(sumClosed, cellI)
    {
        scalar maxOpenness = 0;

        for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
        {
            maxOpenness = max
            (
                maxOpenness,
                mag(sumClosed[cellI][cmpt])
               /(sumMagClosed[cellI][cmpt] + VSMALL)
            );
        }
        openness[cellI] = maxOpenness;

        // Calculate the aspect ration as the maximum of Cartesian component
        // aspect ratio to the total area hydraulic area aspect ratio
        scalar minCmpt = VGREAT;
        scalar maxCmpt = -VGREAT;
        for (direction dir = 0; dir < vector::nComponents; dir++)
        {
            if (meshD[dir] == 1)
            {
                minCmpt = min(minCmpt, sumMagClosed[cellI][dir]);
                maxCmpt = max(maxCmpt, sumMagClosed[cellI][dir]);
            }
        }

        scalar aspectRatio = maxCmpt/(minCmpt + VSMALL);
        if (nDims == 3)
        {
            scalar v = max(VSMALL, vols[cellI]);

            aspectRatio = max
            (
                aspectRatio,
                1.0/6.0*cmptSum(sumMagClosed[cellI])/pow(v, 2.0/3.0)
            );
        }

        aratio[cellI] = aspectRatio;
    }
}


Foam::tmp<Foam::scalarField> Foam::primitiveMeshTools::faceConcavity
(
    const scalar maxSin,
    const primitiveMesh& mesh,
    const pointField& p,
    const vectorField& faceAreas
)
{
    const faceList& fcs = mesh.faces();

    vectorField faceNormals(faceAreas);
    faceNormals /= mag(faceNormals) + VSMALL;

    tmp<scalarField> tfaceAngles(new scalarField(mesh.nFaces()));
    scalarField& faceAngles = tfaceAngles();


    forAll(fcs, faceI)
    {
        const face& f = fcs[faceI];

        // Get edge from f[0] to f[size-1];
        vector ePrev(p[f.first()] - p[f.last()]);
        scalar magEPrev = mag(ePrev);
        ePrev /= magEPrev + VSMALL;

        scalar maxEdgeSin = 0.0;

        forAll(f, fp0)
        {
            // Get vertex after fp
            label fp1 = f.fcIndex(fp0);

            // Normalized vector between two consecutive points
            vector e10(p[f[fp1]] - p[f[fp0]]);
            scalar magE10 = mag(e10);
            e10 /= magE10 + VSMALL;

            if (magEPrev > SMALL && magE10 > SMALL)
            {
                vector edgeNormal = ePrev ^ e10;
                scalar magEdgeNormal = mag(edgeNormal);

                if (magEdgeNormal < maxSin)
                {
                    // Edges (almost) aligned -> face is ok.
                }
                else
                {
                    // Check normal
                    edgeNormal /= magEdgeNormal;

                    if ((edgeNormal & faceNormals[faceI]) < SMALL)
                    {
                        maxEdgeSin = max(maxEdgeSin, magEdgeNormal);
                    }
                }
            }

            ePrev = e10;
            magEPrev = magE10;
        }

        faceAngles[faceI] = maxEdgeSin;
    }

    return tfaceAngles;
}


Foam::tmp<Foam::scalarField> Foam::primitiveMeshTools::faceFlatness
(
    const primitiveMesh& mesh,
    const pointField& p,
    const vectorField& fCtrs,
    const vectorField& faceAreas
)
{
    const faceList& fcs = mesh.faces();

    // Areas are calculated as the sum of areas. (see
    // primitiveMeshFaceCentresAndAreas.C)
    scalarField magAreas(mag(faceAreas));

    tmp<scalarField> tfaceFlatness(new scalarField(mesh.nFaces(), 1.0));
    scalarField& faceFlatness = tfaceFlatness();


    forAll(fcs, faceI)
    {
        const face& f = fcs[faceI];

        if (f.size() > 3 && magAreas[faceI] > VSMALL)
        {
            const point& fc = fCtrs[faceI];

            // Calculate the sum of magnitude of areas and compare to magnitude
            // of sum of areas.

            scalar sumA = 0.0;

            forAll(f, fp)
            {
                const point& thisPoint = p[f[fp]];
                const point& nextPoint = p[f.nextLabel(fp)];

                // Triangle around fc.
                vector n = 0.5*((nextPoint - thisPoint)^(fc - thisPoint));
                sumA += mag(n);
            }

            faceFlatness[faceI] = magAreas[faceI] / (sumA+VSMALL);
        }
    }

    return tfaceFlatness;
}


//Foam::tmp<Foam::scalarField> Foam::primitiveMeshTools::edgeAlignment
//(
//    const primitiveMesh& mesh,
//    const Vector<label>& meshD,
//    const pointField& p
//)
//{
//    label nDirs = 0;
//    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
//    {
//        if (meshD[cmpt] == 1)
//        {
//            nDirs++;
//        }
//        else if (meshD[cmpt] != 0)
//        {
//            FatalErrorIn
//            (
//                "primitiveMeshTools::edgeAlignment"
//                "(const primitiveMesh&, const Vector<label>&"
//                ", const pointField&)"
//            )   << "directions should contain 0 or 1 but is now " << meshD
//                << exit(FatalError);
//        }
//    }
//
//    tmp<scalarField> tedgeAlignment(new scalarField(mesh.nFaces(), 0.0));
//    scalarField& edgeAlignment = tedgeAlignment();
//
//    if (nDirs != vector::nComponents)
//    {
//        const faceList& fcs = mesh.faces();
//
//        forAll(fcs, faceI)
//        {
//            const face& f = fcs[faceI];
//
//            forAll(f, fp)
//            {
//                label p0 = f[fp];
//                label p1 = f.nextLabel(fp);
//                if (p0 < p1)
//                {
//                    vector d(p[p1]-p[p0]);
//                    scalar magD = mag(d);
//
//                    if (magD > ROOTVSMALL)
//                    {
//                        d /= magD;
//
//                        // Check how many empty directions are used by the
//                        // edge.
//                        label nEmptyDirs = 0;
//                        label nNonEmptyDirs = 0;
//                        for
//                        (
//                              direction cmpt=0;
//                              cmpt<vector::nComponents;
//                              cmpt++
//                        )
//                        {
//                            if (mag(d[cmpt]) > 1e-6)
//                            {
//                                if (meshD[cmpt] == 0)
//                                {
//                                    nEmptyDirs++;
//                                }
//                                else
//                                {
//                                    nNonEmptyDirs++;
//                                }
//                            }
//                        }
//
//                        if (nEmptyDirs == 0)
//                        {
//                            // Purely in ok directions.
//                        }
//                        else if (nEmptyDirs == 1)
//                        {
//                            // Ok if purely in empty directions.
//                            if (nNonEmptyDirs > 0)
//                            {
//                                edgeAlignment[faceI] = max
//                                (
//                                    edgeAlignment[faceI],
//                                    1
//                                );
//                            }
//                        }
//                        else if (nEmptyDirs > 1)
//                        {
//                            // Always an error
//                            edgeAlignment[faceI] = max
//                            (
//                                edgeAlignment[faceI],
//                                2
//                            );
//                        }
//                    }
//                }
//            }
//        }
//    }
//
//    return tedgeAlignment;
//}


// Checks cells with 1 or less internal faces. Give numerical problems.
Foam::tmp<Foam::scalarField> Foam::primitiveMeshTools::cellDeterminant
(
    const primitiveMesh& mesh,
    const Vector<label>& meshD,
    const vectorField& faceAreas,
    const PackedBoolList& internalOrCoupledFace
)
{
    // Determine number of dimensions and (for 2D) missing dimension
    label nDims = 0;
    label twoD = -1;
    for (direction dir = 0; dir < vector::nComponents; dir++)
    {
        if (meshD[dir] == 1)
        {
            nDims++;
        }
        else
        {
            twoD = dir;
        }
    }



    tmp<scalarField> tcellDeterminant(new scalarField(mesh.nCells()));
    scalarField& cellDeterminant = tcellDeterminant();

    const cellList& c = mesh.cells();

    if (nDims == 1)
    {
        cellDeterminant = 1.0;
    }
    else
    {
        forAll (c, cellI)
        {
            const labelList& curFaces = c[cellI];

            // Calculate local normalization factor
            scalar avgArea = 0;

            label nInternalFaces = 0;

            forAll(curFaces, i)
            {
                if (internalOrCoupledFace[curFaces[i]])
                {
                    avgArea += mag(faceAreas[curFaces[i]]);

                    nInternalFaces++;
                }
            }

            if (nInternalFaces == 0)
            {
                cellDeterminant[cellI] = 0;
            }
            else
            {
                avgArea /= nInternalFaces;

                symmTensor areaTensor(symmTensor::zero);

                forAll(curFaces, i)
                {
                    if (internalOrCoupledFace[curFaces[i]])
                    {
                        areaTensor += sqr(faceAreas[curFaces[i]]/avgArea);
                    }
                }

                if (nDims == 2)
                {
                    // Add the missing eigenvector (such that it does not
                    // affect the determinant)
                    if (twoD == 0)
                    {
                        areaTensor.xx() = 1;
                    }
                    else if (twoD == 1)
                    {
                        areaTensor.yy() = 1;
                    }
                    else
                    {
                        areaTensor.zz() = 1;
                    }
                }

                cellDeterminant[cellI] = mag(det(areaTensor));
            }
        }
    }

    return tcellDeterminant;
}


// ************************************************************************* //
