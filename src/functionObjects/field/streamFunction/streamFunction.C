/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2019 OpenFOAM Foundation
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

#include "streamFunction.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "emptyPolyPatch.H"
#include "symmetryPlanePolyPatch.H"
#include "symmetryPolyPatch.H"
#include "wedgePolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(streamFunction, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        streamFunction,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::pointScalarField> Foam::functionObjects::streamFunction::calc
(
    const surfaceScalarField& phi
) const
{
    Log << "    functionObjects::" << type() << " " << name()
        << " calculating steam-function" << endl;

    Vector<label> slabNormal((Vector<label>::one - mesh_.geometricD())/2);
    const direction slabDir
    (
        slabNormal
      & Vector<label>(Vector<label>::X, Vector<label>::Y, Vector<label>::Z)
    );

    scalar thickness = vector(slabNormal) & mesh_.bounds().span();

    const pointMesh& pMesh = pointMesh::New(mesh_);

    tmp<pointScalarField> tstreamFunction
    (
        pointScalarField::New
        (
            "streamFunction",
            pMesh,
            dimensionedScalar(phi.dimensions(), 0)
        )
    );
    pointScalarField& streamFunction = tstreamFunction.ref();

    labelList visitedPoint(mesh_.nPoints());
    forAll(visitedPoint, pointi)
    {
        visitedPoint[pointi] = 0;
    }
    label nVisited = 0;
    label nVisitedOld = 0;

    const faceUList& faces = mesh_.faces();
    const pointField& points = mesh_.points();

    label nInternalFaces = mesh_.nInternalFaces();

    vectorField unitAreas(mesh_.faceAreas());
    unitAreas /= mag(unitAreas);

    const polyPatchList& patches = mesh_.boundaryMesh();

    bool finished = true;

    // Find the boundary face with zero flux. set the stream function
    // to zero on that face
    bool found = false;

    do
    {
        found = false;

        forAll(patches, patchi)
        {
            const primitivePatch& bouFaces = patches[patchi];

            if (!isType<emptyPolyPatch>(patches[patchi]))
            {
                forAll(bouFaces, facei)
                {
                    if
                    (
                        magSqr(phi.boundaryField()[patchi][facei]) < small
                    )
                    {
                        const labelList& zeroPoints = bouFaces[facei];

                        // Zero flux face found
                        found = true;

                        forAll(zeroPoints, pointi)
                        {
                            if (visitedPoint[zeroPoints[pointi]] == 1)
                            {
                                found = false;
                                break;
                            }
                        }

                        if (found)
                        {
                            Log << "        Zero face: patch: " << patchi
                                << " face: " << facei << endl;

                            forAll(zeroPoints, pointi)
                            {
                                streamFunction[zeroPoints[pointi]] = 0;
                                visitedPoint[zeroPoints[pointi]] = 1;
                                nVisited++;
                            }

                            break;
                        }
                    }
                }
            }

            if (found) break;
        }

        if (!found)
        {
            Log << "        Zero flux boundary face not found. "
                << "Using cell as a reference."
                << endl;

            const cellList& c = mesh_.cells();

            forAll(c, cI)
            {
                labelList zeroPoints = c[cI].labels(mesh_.faces());

                bool found = true;

                forAll(zeroPoints, pointi)
                {
                    if (visitedPoint[zeroPoints[pointi]] == 1)
                    {
                        found = false;
                        break;
                    }
                }

                if (found)
                {
                    forAll(zeroPoints, pointi)
                    {
                        streamFunction[zeroPoints[pointi]] = 0.0;
                        visitedPoint[zeroPoints[pointi]] = 1;
                        nVisited++;
                    }

                    break;
                }
                else
                {
                    FatalErrorInFunction
                        << "Cannot find initialisation face or a cell."
                        << exit(FatalError);
                }
            }
        }

        // Loop through all faces. If one of the points on
        // the face has the streamfunction value different
        // from -1, all points with -1 ont that face have the
        // streamfunction value equal to the face flux in
        // that point plus the value in the visited point
        do
        {
            finished = true;

            for (label facei = nInternalFaces; facei<faces.size(); facei++)
            {
                const labelList& curBPoints = faces[facei];
                bool bPointFound = false;

                scalar currentBStream = 0.0;
                vector currentBStreamPoint(0, 0, 0);

                forAll(curBPoints, pointi)
                {
                    // Check if the point has been visited
                    if (visitedPoint[curBPoints[pointi]] == 1)
                    {
                        // The point has been visited
                        currentBStream = streamFunction[curBPoints[pointi]];
                        currentBStreamPoint = points[curBPoints[pointi]];

                        bPointFound = true;

                        break;
                    }
                }

                if (bPointFound)
                {
                    // Sort out other points on the face
                    forAll(curBPoints, pointi)
                    {
                        // Check if the point has been visited
                        if (visitedPoint[curBPoints[pointi]] == 0)
                        {
                            label patchNo =
                                mesh_.boundaryMesh().whichPatch(facei);

                            if
                            (
                                !isType<emptyPolyPatch>(patches[patchNo])
                             && !isType<symmetryPlanePolyPatch>
                                 (patches[patchNo])
                             && !isType<symmetryPolyPatch>(patches[patchNo])
                             && !isType<wedgePolyPatch>(patches[patchNo])
                            )
                            {
                                label faceNo =
                                    mesh_.boundaryMesh()[patchNo]
                                   .whichFace(facei);

                                vector edgeHat =
                                    points[curBPoints[pointi]]
                                  - currentBStreamPoint;
                                edgeHat.replace(slabDir, 0);
                                edgeHat /= mag(edgeHat);

                                vector nHat = unitAreas[facei];

                                if (edgeHat.y() > vSmall)
                                {
                                    visitedPoint[curBPoints[pointi]] = 1;
                                    nVisited++;

                                    streamFunction[curBPoints[pointi]] =
                                        currentBStream
                                      + phi.boundaryField()[patchNo][faceNo]
                                       *sign(nHat.x());
                                }
                                else if (edgeHat.y() < -vSmall)
                                {
                                    visitedPoint[curBPoints[pointi]] = 1;
                                    nVisited++;

                                    streamFunction[curBPoints[pointi]] =
                                        currentBStream
                                      - phi.boundaryField()[patchNo][faceNo]
                                       *sign(nHat.x());
                                }
                                else
                                {
                                    if (edgeHat.x() > vSmall)
                                    {
                                        visitedPoint[curBPoints[pointi]] = 1;
                                        nVisited++;

                                        streamFunction[curBPoints[pointi]] =
                                            currentBStream
                                          + phi.boundaryField()[patchNo][faceNo]
                                           *sign(nHat.y());
                                    }
                                    else if (edgeHat.x() < -vSmall)
                                    {
                                        visitedPoint[curBPoints[pointi]] = 1;
                                        nVisited++;

                                        streamFunction[curBPoints[pointi]] =
                                            currentBStream
                                          - phi.boundaryField()[patchNo][faceNo]
                                           *sign(nHat.y());
                                    }
                                }
                            }
                        }
                    }
                }
                else
                {
                    finished = false;
                }
            }

            for (label facei=0; facei<nInternalFaces; facei++)
            {
                // Get the list of point labels for the face
                const labelList& curPoints = faces[facei];

                bool pointFound = false;

                scalar currentStream = 0.0;
                point currentStreamPoint(0, 0, 0);

                forAll(curPoints, pointi)
                {
                    // Check if the point has been visited
                    if (visitedPoint[curPoints[pointi]] == 1)
                    {
                        // The point has been visited
                        currentStream = streamFunction[curPoints[pointi]];
                        currentStreamPoint = points[curPoints[pointi]];
                        pointFound = true;

                        break;
                    }
                }

                if (pointFound)
                {
                    // Sort out other points on the face
                    forAll(curPoints, pointi)
                    {
                        // Check if the point has been visited
                        if (visitedPoint[curPoints[pointi]] == 0)
                        {
                            vector edgeHat =
                                points[curPoints[pointi]] - currentStreamPoint;

                            edgeHat.replace(slabDir, 0);
                            edgeHat /= mag(edgeHat);

                            vector nHat = unitAreas[facei];

                            if (edgeHat.y() > vSmall)
                            {
                                visitedPoint[curPoints[pointi]] = 1;
                                nVisited++;

                                streamFunction[curPoints[pointi]] =
                                  currentStream
                                + phi[facei]*sign(nHat.x());
                            }
                            else if (edgeHat.y() < -vSmall)
                            {
                                visitedPoint[curPoints[pointi]] = 1;
                                nVisited++;

                                streamFunction[curPoints[pointi]] =
                                  currentStream
                                - phi[facei]*sign(nHat.x());
                            }
                        }
                    }
                }
                else
                {
                    finished = false;
                }
            }

            if (nVisited == nVisitedOld)
            {
                // Find new seed.  This must be a
                // multiply connected domain
                Log << "        Exhausted a seed, looking for new seed "
                    << "(this is correct for multiply connected domains).";

                break;
            }
            else
            {
                nVisitedOld = nVisited;
            }
        } while (!finished);
    } while (!finished);

    // Normalise the stream-function by the 2D mesh thickness
    streamFunction /= thickness;
    streamFunction.boundaryFieldRef() = 0.0;

    return tstreamFunction;
}


bool Foam::functionObjects::streamFunction::calc()
{
    if (foundObject<surfaceScalarField>(fieldName_))
    {
        const surfaceScalarField& phi =
            mesh_.lookupObject<surfaceScalarField>(fieldName_);

        return store(resultName_, calc(phi));
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::streamFunction::streamFunction
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, typeName, "phi")
{
    label nD = mesh_.nGeometricD();

    if (nD != 2)
    {
        FatalErrorInFunction
            << "Case is not 2D, stream-function cannot be computed"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::streamFunction::~streamFunction()
{}


// ************************************************************************* //
