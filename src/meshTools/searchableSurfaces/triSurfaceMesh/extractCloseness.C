/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2021 OpenFOAM Foundation
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

#include "triSurfaceMesh.H"
#include "triSurfaceFields.H"
#include "meshTools.H"
#include "Time.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::triSurfaceMesh::drawHitProblem
(
    const label fi,
    const point& start,
    const point& p,
    const point& end,
    const pointIndexHitList& hitInfo
) const
{
    if (debug)
    {
        const List<labelledTri>& tris = *this;
        const pointField& points = this->points();

        Info<< nl << "# findLineAll did not hit its own face."
            << nl << "# fi " << fi
            << nl << "# start " << start
            << nl << "# point " << p
            << nl << "# end " << end
            << nl << "# hitInfo " << hitInfo
            << endl;

        meshTools::writeOBJ(Info, start);
        meshTools::writeOBJ(Info, p);
        meshTools::writeOBJ(Info, end);

        Info<< "l 1 2 3" << endl;

        meshTools::writeOBJ(Info, points[tris[fi][0]]);
        meshTools::writeOBJ(Info, points[tris[fi][1]]);
        meshTools::writeOBJ(Info, points[tris[fi][2]]);

        Info<< "f 4 5 6" << endl;

        forAll(hitInfo, hi)
        {
            label hfi = hitInfo[hi].index();

            meshTools::writeOBJ(Info, points[tris[hfi][0]]);
            meshTools::writeOBJ(Info, points[tris[hfi][1]]);
            meshTools::writeOBJ(Info, points[tris[hfi][2]]);

            Info<< "f "
                << 3*hi + 7 << " "
                << 3*hi + 8 << " "
                << 3*hi + 9
                << endl;
        }
    }
}


void Foam::triSurfaceMesh::processHit
(
    scalar& internalCloseness,
    scalar& externalCloseness,
    const scalar internalToleranceCosAngle,
    const scalar externalToleranceCosAngle,
    const label fi,
    const point& start,
    const point& p,
    const point& end,
    const vector& normal,
    const vectorField& normals,
    const pointIndexHitList& hitInfo
) const
{
    if (hitInfo.size() < 1)
    {
        drawHitProblem(fi, start, p, end, hitInfo);
    }
    else if (hitInfo.size() == 1)
    {
        if (!hitInfo[0].hit())
        {
        }
        else if (hitInfo[0].index() != fi)
        {
            drawHitProblem(fi, start, p, end, hitInfo);
        }
    }
    else
    {
        label ownHiti = -1;

        forAll(hitInfo, hI)
        {
            // Find the hit on the triangle that launched the ray

            if (hitInfo[hI].index() == fi)
            {
                ownHiti = hI;
                break;
            }
        }

        if (ownHiti < 0)
        {
            drawHitProblem(fi, start, p, end, hitInfo);
        }
        else if (ownHiti == 0)
        {
            // There are no internal hits, the first hit is the
            // closest external hit

            if
            (
                (normal & normals[hitInfo[ownHiti + 1].index()])
              < externalToleranceCosAngle
            )
            {
                externalCloseness = min
                (
                    externalCloseness,
                    mag(p - hitInfo[ownHiti + 1].hitPoint())
                );
            }
        }
        else if (ownHiti == hitInfo.size() - 1)
        {
            // There are no external hits, the last but one hit is
            // the closest internal hit

            if
            (
                (normal & normals[hitInfo[ownHiti - 1].index()])
              < internalToleranceCosAngle
            )
            {
                internalCloseness = min
                (
                    internalCloseness,
                    mag(p - hitInfo[ownHiti - 1].hitPoint())
                );
            }
        }
        else
        {
            if
            (
                (normal & normals[hitInfo[ownHiti + 1].index()])
              < externalToleranceCosAngle
            )
            {
                externalCloseness = min
                (
                    externalCloseness,
                    mag(p - hitInfo[ownHiti + 1].hitPoint())
                );
            }

            if
            (
                (normal & normals[hitInfo[ownHiti - 1].index()])
              < internalToleranceCosAngle
            )
            {
                internalCloseness = min
                (
                    internalCloseness,
                    mag(p - hitInfo[ownHiti - 1].hitPoint())
                );
            }
        }
    }
}


Foam::Pair<Foam::tmp<Foam::triSurfaceScalarField>>
Foam::triSurfaceMesh::extractCloseness
(
    const scalar internalAngleTolerance,
    const scalar externalAngleTolerance
) const
{
    const scalar internalToleranceCosAngle
    (
        cos(degToRad(180 - internalAngleTolerance))
    );

    const scalar externalToleranceCosAngle
    (
        cos(degToRad(180 - externalAngleTolerance))
    );

    const Time& runTime = objectRegistry::time();

    // Prepare start and end points for intersection tests

    const vectorField& normals = faceNormals();

    const scalar span = bounds().mag();

    const pointField start(faceCentres() - span*normals);
    const pointField end(faceCentres() + span*normals);
    const pointField& faceCentres = this->faceCentres();

    List<List<pointIndexHit>> allHitinfo;

    // Find all intersections (in order)
    findLineAll(start, end, allHitinfo);

    scalarField internalCloseness(start.size(), great);
    scalarField externalCloseness(start.size(), great);

    forAll(allHitinfo, fi)
    {
        const List<pointIndexHit>& hitInfo = allHitinfo[fi];

        processHit
        (
            internalCloseness[fi],
            externalCloseness[fi],
            internalToleranceCosAngle,
            externalToleranceCosAngle,
            fi,
            start[fi],
            faceCentres[fi],
            end[fi],
            normals[fi],
            normals,
            hitInfo
        );
    }

    return Pair<tmp<triSurfaceScalarField>>
    (
        tmp<triSurfaceScalarField>
        (
            new triSurfaceScalarField
            (
                IOobject
                (
                    objectRegistry::name() + ".internalCloseness",
                    runTime.constant(),
                    searchableSurface::geometryDir(runTime),
                    runTime
                ),
                *this,
                dimLength,
                internalCloseness
            )
        ),

        tmp<triSurfaceScalarField>
        (
            new triSurfaceScalarField
            (
                IOobject
                (
                    objectRegistry::name() + ".externalCloseness",
                    runTime.constant(),
                    searchableSurface::geometryDir(runTime),
                    runTime
                ),
                *this,
                dimLength,
                externalCloseness
            )
        )
    );
}


Foam::Pair<Foam::tmp<Foam::triSurfacePointScalarField>>
Foam::triSurfaceMesh::extractPointCloseness
(
    const scalar internalAngleTolerance,
    const scalar externalAngleTolerance
) const
{
    const scalar internalToleranceCosAngle
    (
        cos(degToRad(180 - internalAngleTolerance))
    );

    const scalar externalToleranceCosAngle
    (
        cos(degToRad(180 - externalAngleTolerance))
    );

    const Time& runTime = objectRegistry::time();

    // Prepare start and end points for intersection tests

    const pointField& points = this->points();
    const labelList& meshPoints = this->meshPoints();
    const pointField& faceCentres = this->faceCentres();
    const vectorField& normals = this->faceNormals();
    const labelListList& pointFaces = this->pointFaces();

    const scalar span = bounds().mag();

    label nPointFaces = 0;
    forAll(pointFaces, pfi)
    {
        nPointFaces += pointFaces[pfi].size();
    }

    pointField facePoints(nPointFaces);
    pointField start(nPointFaces);
    pointField end(nPointFaces);

    label i = 0;
    forAll(points, pi)
    {
        forAll(pointFaces[pi], pfi)
        {
            const label fi = pointFaces[pi][pfi];

            facePoints[i] = (0.9*points[meshPoints[pi]] + 0.1*faceCentres[fi]);
            const vector& n = normals[fi];

            start[i] = facePoints[i] - span*n;
            end[i] = facePoints[i] + span*n;

            i++;
        }
    }

    List<pointIndexHitList> allHitinfo;

    // Find all intersections (in order)
    findLineAll(start, end, allHitinfo);

    scalarField internalCloseness(points.size(), great);
    scalarField externalCloseness(points.size(), great);

    i = 0;
    forAll(points, pi)
    {
        forAll(pointFaces[pi], pfi)
        {
            const label fi = pointFaces[pi][pfi];
            const pointIndexHitList& hitInfo = allHitinfo[i];

            processHit
            (
                internalCloseness[pi],
                externalCloseness[pi],
                internalToleranceCosAngle,
                externalToleranceCosAngle,
                fi,
                start[i],
                facePoints[i],
                end[i],
                normals[fi],
                normals,
                hitInfo
            );

            i++;
        }
    }

    return Pair<tmp<triSurfacePointScalarField>>
    (
        tmp<triSurfacePointScalarField>
        (
            new triSurfacePointScalarField
            (
                IOobject
                (
                    objectRegistry::name() + ".internalPointCloseness",
                    runTime.constant(),
                    searchableSurface::geometryDir(runTime),
                    runTime
                ),
                *this,
                dimLength,
                internalCloseness
            )
        ),

        tmp<triSurfacePointScalarField>
        (
            new triSurfacePointScalarField
            (
                IOobject
                (
                    objectRegistry::name() + ".externalPointCloseness",
                    runTime.constant(),
                    searchableSurface::geometryDir(runTime),
                    runTime
                ),
                *this,
                dimLength,
                externalCloseness
            )
        )
    );
}


// ************************************************************************* //
