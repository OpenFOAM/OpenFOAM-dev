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

\*---------------------------------------------------------------------------*/

#include "motionSmootherAlgo.H"
#include "polyMeshCheck.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::motionSmootherAlgo::checkMesh
(
    const bool report,
    const polyMesh& mesh,
    const dictionary& dict,
    const labelList& checkFaces,
    labelHashSet& wrongFaces
)
{
    List<labelPair> emptyBaffles;
    return checkMesh
    (
        report,
        mesh,
        dict,
        checkFaces,
        emptyBaffles,
        wrongFaces
    );
}

bool Foam::motionSmootherAlgo::checkMesh
(
    const bool report,
    const polyMesh& mesh,
    const dictionary& dict,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet& wrongFaces
)
{
    const scalar maxNonOrtho
    (
        dict.lookup<scalar>("maxNonOrtho", true)
    );
    const scalar minVol
    (
        dict.lookup<scalar>("minVol", true)
    );
    const scalar minTetQuality
    (
        dict.lookup<scalar>("minTetQuality", true)
    );
    const scalar maxConcave
    (
        dict.lookup<scalar>("maxConcave", true)
    );
    const scalar maxIntSkew
    (
        dict.lookup<scalar>("maxInternalSkewness", true)
    );
    const scalar maxBounSkew
    (
        dict.lookup<scalar>("maxBoundarySkewness", true)
    );
    const scalar minWeight
    (
        dict.lookup<scalar>("minFaceWeight", true)
    );
    const scalar minVolRatio
    (
        dict.lookup<scalar>("minVolRatio", true)
    );
    const scalar minTwist
    (
        dict.lookup<scalar>("minTwist", true)
    );
    scalar minFaceFlatness = -1.0;
    dict.readIfPresent("minFaceFlatness", minFaceFlatness, true);
    const scalar minDet
    (
        dict.lookup<scalar>("minDeterminant", true)
    );
    label nWrongFaces = 0;

    Info<< "Checking faces in error :" << endl;
    // Pout.setf(ios_base::left);

    if (maxNonOrtho < 180.0-small)
    {
        polyMeshCheck::checkFaceDotProduct
        (
            report,
            maxNonOrtho,
            mesh,
            mesh.cellCentres(),
            mesh.faceAreas(),
            checkFaces,
            baffles,
            &wrongFaces
        );

        label nNewWrongFaces = returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    non-orthogonality > "
            << setw(3) << maxNonOrtho
            << " degrees                        : "
            << nNewWrongFaces-nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (minVol > -great)
    {
        const scalar refVol = pow3(mesh.bounds().minDim());

        polyMeshCheck::checkFacePyramids
        (
            report,
            minVol*refVol,
            mesh,
            mesh.cellCentres(),
            mesh.points(),
            checkFaces,
            baffles,
            &wrongFaces
        );

        label nNewWrongFaces = returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    faces with face pyramid volume < "
            << setw(5) << minVol << "                 : "
            << nNewWrongFaces-nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (minTetQuality > -great)
    {
        polyMeshCheck::checkFaceTets
        (
            report,
            minTetQuality,
            mesh,
            mesh.cellCentres(),
            mesh.faceCentres(),
            mesh.points(),
            checkFaces,
            baffles,
            &wrongFaces
        );

        label nNewWrongFaces = returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    faces with face-decomposition tet quality < "
            << setw(5) << minTetQuality << "      : "
            << nNewWrongFaces-nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (maxConcave < 180.0-small)
    {
        polyMeshCheck::checkFaceAngles
        (
            report,
            maxConcave,
            mesh,
            mesh.faceAreas(),
            mesh.points(),
            checkFaces,
            &wrongFaces
        );

        label nNewWrongFaces = returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    faces with concavity > "
            << setw(3) << maxConcave
            << " degrees                     : "
            << nNewWrongFaces-nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (maxIntSkew > 0 || maxBounSkew > 0)
    {
        polyMeshCheck::checkFaceSkewness
        (
            report,
            maxIntSkew,
            maxBounSkew,
            mesh,
            mesh.points(),
            mesh.cellCentres(),
            mesh.faceCentres(),
            mesh.faceAreas(),
            checkFaces,
            baffles,
            &wrongFaces
        );

        label nNewWrongFaces = returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    faces with skewness > "
            << setw(3) << maxIntSkew
            << " (internal) or " << setw(3) << maxBounSkew
            << " (boundary) : " << nNewWrongFaces-nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (minWeight >= 0 && minWeight < 1)
    {
        polyMeshCheck::checkFaceWeights
        (
            report,
            minWeight,
            mesh,
            mesh.cellCentres(),
            mesh.faceCentres(),
            mesh.faceAreas(),
            checkFaces,
            baffles,
            &wrongFaces
        );

        label nNewWrongFaces = returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    faces with interpolation weights (0..1)  < "
            << setw(5) << minWeight
            << "       : "
            << nNewWrongFaces-nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (minVolRatio >= 0)
    {
        polyMeshCheck::checkVolRatio
        (
            report,
            minVolRatio,
            mesh,
            mesh.cellVolumes(),
            checkFaces,
            baffles,
            &wrongFaces
        );

        label nNewWrongFaces = returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    faces with volume ratio of neighbour cells < "
            << setw(5) << minVolRatio
            << "     : "
            << nNewWrongFaces-nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (minTwist > -1)
    {
        // Pout<< "Checking face twist: dot product of face normal "
        //    << "with face triangle normals" << endl;
        polyMeshCheck::checkFaceTwist
        (
            report,
            minTwist,
            mesh,
            mesh.cellCentres(),
            mesh.faceAreas(),
            mesh.faceCentres(),
            mesh.points(),
            checkFaces,
            &wrongFaces
        );

        label nNewWrongFaces = returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    faces with face twist < "
            << setw(5) << minTwist
            << "                          : "
            << nNewWrongFaces-nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (minFaceFlatness > -small)
    {
        polyMeshCheck::checkFaceFlatness
        (
            report,
            minFaceFlatness,
            mesh,
            mesh.faceAreas(),
            mesh.faceCentres(),
            mesh.points(),
            checkFaces,
            &wrongFaces
        );

        label nNewWrongFaces = returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    faces with flatness < "
            << setw(5) << minFaceFlatness
            << "                      : "
            << nNewWrongFaces-nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    if (minDet > -1)
    {
        polyMeshCheck::checkCellDeterminant
        (
            report,
            minDet,
            mesh,
            mesh.faceAreas(),
            checkFaces,
            &wrongFaces
        );

        label nNewWrongFaces = returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    faces on cells with determinant < "
            << setw(5) << minDet << "                : "
            << nNewWrongFaces-nWrongFaces << endl;

        nWrongFaces = nNewWrongFaces;
    }

    // Pout.setf(ios_base::right);

    return nWrongFaces > 0;
}


bool Foam::motionSmootherAlgo::checkMesh
(
    const bool report,
    const polyMesh& mesh,
    const dictionary& dict,
    labelHashSet& wrongFaces
)
{
    return checkMesh
    (
        report,
        mesh,
        dict,
        identityMap(mesh.nFaces()),
        wrongFaces
    );
}


// ************************************************************************* //
