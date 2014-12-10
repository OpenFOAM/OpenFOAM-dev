/*---------------------------------------------------------------------------* \
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

Application
    streamFunction

Description
    Calculates and writes the stream function of velocity field U at each
    time.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pointFields.H"
#include "emptyPolyPatch.H"
#include "symmetryPlanePolyPatch.H"
#include "symmetryPolyPatch.H"
#include "wedgePolyPatch.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "addRegionOption.H"

    #include "setRootCase.H"
    #include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    #include "createNamedMesh.H"

    label nD = mesh.nGeometricD();

    if (nD != 2)
    {
        FatalErrorIn(args.executable())
            << "Case is not 2D, stream-function cannot be computed"
            << exit(FatalError);
    }

    const vector slabDir((Vector<label>::one - mesh.geometricD())/2);
    scalar thickness = slabDir & mesh.bounds().span();

    const pointMesh& pMesh = pointMesh::New(mesh);

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< nl << "Time: " << runTime.timeName() << endl;

        IOobject phiHeader
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ
        );

        if (phiHeader.headerOk())
        {
            mesh.readUpdate();

            Info<< nl << "Reading field phi" << endl;

            surfaceScalarField phi
            (
                IOobject
                (
                    "phi",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );

            pointScalarField streamFunction
            (
                IOobject
                (
                    "streamFunction",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                pMesh,
                dimensionedScalar("zero", phi.dimensions(), 0.0)
            );

            labelList visitedPoint(mesh.nPoints());
            forAll(visitedPoint, pointI)
            {
                visitedPoint[pointI] = 0;
            }
            label nVisited = 0;
            label nVisitedOld = 0;

            const faceUList& faces = mesh.faces();
            const pointField& points = mesh.points();

            label nInternalFaces = mesh.nInternalFaces();

            vectorField unitAreas(mesh.faceAreas());
            unitAreas /= mag(unitAreas);

            const polyPatchList& patches = mesh.boundaryMesh();

            bool finished = true;

            // Find the boundary face with zero flux. set the stream function
            // to zero on that face
            bool found = false;

            do
            {
                found = false;

                forAll(patches, patchI)
                {
                    const primitivePatch& bouFaces = patches[patchI];

                    if (!isType<emptyPolyPatch>(patches[patchI]))
                    {
                        forAll(bouFaces, faceI)
                        {
                            if
                            (
                                magSqr(phi.boundaryField()[patchI][faceI])
                              < SMALL
                            )
                            {
                                const labelList& zeroPoints = bouFaces[faceI];

                                // Zero flux face found
                                found = true;

                                forAll(zeroPoints, pointI)
                                {
                                    if (visitedPoint[zeroPoints[pointI]] == 1)
                                    {
                                        found = false;
                                        break;
                                    }
                                }

                                if (found)
                                {
                                    Info<< "Zero face: patch: " << patchI
                                        << "    face: " << faceI << endl;

                                    forAll(zeroPoints, pointI)
                                    {
                                        streamFunction[zeroPoints[pointI]] = 0;
                                        visitedPoint[zeroPoints[pointI]] = 1;
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
                    Info<< "zero flux boundary face not found. "
                        << "Using cell as a reference."
                        << endl;

                    const cellList& c = mesh.cells();

                    forAll(c, cI)
                    {
                        labelList zeroPoints = c[cI].labels(mesh.faces());

                        bool found = true;

                        forAll(zeroPoints, pointI)
                        {
                            if (visitedPoint[zeroPoints[pointI]] == 1)
                            {
                                found = false;
                                break;
                            }
                        }

                        if (found)
                        {
                            forAll(zeroPoints, pointI)
                            {
                                streamFunction[zeroPoints[pointI]] = 0.0;
                                visitedPoint[zeroPoints[pointI]] = 1;
                                nVisited++;
                            }

                            break;
                        }
                        else
                        {
                            FatalErrorIn(args.executable())
                                << "Cannot find initialisation face or a cell."
                                << abort(FatalError);
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

                     for
                     (
                         label faceI = nInternalFaces;
                         faceI<faces.size();
                         faceI++
                     )
                     {
                         const labelList& curBPoints = faces[faceI];
                         bool bPointFound = false;

                         scalar currentBStream = 0.0;
                         vector currentBStreamPoint(0, 0, 0);

                         forAll(curBPoints, pointI)
                         {
                             // Check if the point has been visited
                             if (visitedPoint[curBPoints[pointI]] == 1)
                             {
                                 // The point has been visited
                                 currentBStream =
                                     streamFunction[curBPoints[pointI]];
                                 currentBStreamPoint =
                                     points[curBPoints[pointI]];

                                 bPointFound = true;

                                 break;
                             }
                         }

                         if (bPointFound)
                         {
                             // Sort out other points on the face
                             forAll(curBPoints, pointI)
                             {
                                 // Check if the point has been visited
                                 if (visitedPoint[curBPoints[pointI]] == 0)
                                 {
                                     label patchNo =
                                         mesh.boundaryMesh().whichPatch(faceI);

                                     if
                                     (
                                        !isType<emptyPolyPatch>
                                         (patches[patchNo])
                                     && !isType<symmetryPlanePolyPatch>
                                         (patches[patchNo])
                                     && !isType<symmetryPolyPatch>
                                         (patches[patchNo])
                                     && !isType<wedgePolyPatch>
                                         (patches[patchNo])
                                     )
                                     {
                                         label faceNo =
                                             mesh.boundaryMesh()[patchNo]
                                             .whichFace(faceI);

                                         vector edgeHat =
                                             points[curBPoints[pointI]]
                                             - currentBStreamPoint;
                                         edgeHat.replace(vector::Z, 0);
                                         edgeHat /= mag(edgeHat);

                                         vector nHat = unitAreas[faceI];

                                         if (edgeHat.y() > VSMALL)
                                         {
                                             visitedPoint[curBPoints[pointI]] =
                                                 1;
                                             nVisited++;

                                             streamFunction[curBPoints[pointI]]
                                                 =
                                                 currentBStream
                                               + phi.boundaryField()
                                                 [patchNo][faceNo]
                                                 *sign(nHat.x());
                                         }
                                         else if (edgeHat.y() < -VSMALL)
                                         {
                                             visitedPoint[curBPoints[pointI]] =
                                                 1;
                                             nVisited++;

                                             streamFunction[curBPoints[pointI]]
                                                 =
                                                 currentBStream
                                               - phi.boundaryField()
                                                 [patchNo][faceNo]
                                                 *sign(nHat.x());
                                         }
                                         else
                                         {
                                             if (edgeHat.x() > VSMALL)
                                             {
                                                 visitedPoint
                                                     [curBPoints[pointI]] = 1;
                                                 nVisited++;

                                                 streamFunction
                                                     [curBPoints[pointI]] =
                                                     currentBStream
                                                   + phi.boundaryField()
                                                     [patchNo][faceNo]
                                                     *sign(nHat.y());
                                             }
                                             else if (edgeHat.x() < -VSMALL)
                                             {
                                                 visitedPoint
                                                     [curBPoints[pointI]] = 1;
                                                 nVisited++;

                                                 streamFunction
                                                     [curBPoints[pointI]] =
                                                     currentBStream
                                                   - phi.boundaryField()
                                                     [patchNo][faceNo]
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

                     for (label faceI=0; faceI<nInternalFaces; faceI++)
                     {
                         // Get the list of point labels for the face
                         const labelList& curPoints = faces[faceI];

                         bool pointFound = false;

                         scalar currentStream = 0.0;
                         point currentStreamPoint(0, 0, 0);

                         forAll(curPoints, pointI)
                         {
                             // Check if the point has been visited
                             if (visitedPoint[curPoints[pointI]] == 1)
                             {
                                 // The point has been visited
                                 currentStream =
                                     streamFunction[curPoints[pointI]];
                                 currentStreamPoint =
                                     points[curPoints[pointI]];
                                 pointFound = true;

                                 break;
                             }
                         }

                         if (pointFound)
                         {
                             // Sort out other points on the face
                             forAll(curPoints, pointI)
                             {
                                 // Check if the point has been visited
                                 if (visitedPoint[curPoints[pointI]] == 0)
                                 {
                                     vector edgeHat =
                                         points[curPoints[pointI]]
                                       - currentStreamPoint;

                                     edgeHat.replace(vector::Z, 0);
                                     edgeHat /= mag(edgeHat);

                                     vector nHat = unitAreas[faceI];

                                     if (edgeHat.y() > VSMALL)
                                     {
                                         visitedPoint[curPoints[pointI]] = 1;
                                         nVisited++;

                                         streamFunction[curPoints[pointI]] =
                                             currentStream
                                           + phi[faceI]*sign(nHat.x());
                                     }
                                     else if (edgeHat.y() < -VSMALL)
                                     {
                                         visitedPoint[curPoints[pointI]] = 1;
                                         nVisited++;

                                         streamFunction[curPoints[pointI]] =
                                             currentStream
                                           - phi[faceI]*sign(nHat.x());
                                     }
                                 }
                             }
                         }
                         else
                         {
                             finished = false;
                         }
                     }

                     Info<< ".";

                     if (nVisited == nVisitedOld)
                     {
                         // Find new seed.  This must be a
                         // multiply connected domain
                         Info<< nl << "Exhausted a seed. Looking for new seed "
                             << "(this is correct for multiply connected "
                             << "domains).";

                         break;
                     }
                     else
                     {
                         nVisitedOld = nVisited;
                     }
                } while (!finished);

                Info<< endl;
            } while (!finished);

            // Normalise the stream-function by the 2D mesh thickness
            streamFunction /= thickness;
            streamFunction.boundaryField() = 0.0;
            streamFunction.write();
        }
        else
        {
            WarningIn(args.executable())
                << "Flux field does not exist."
                << " Stream function not calculated" << endl;
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
