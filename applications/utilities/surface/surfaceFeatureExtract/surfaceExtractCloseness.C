/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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

#include "surfaceFeatureExtract.H"
#include "Time.H"
#include "triSurfaceMesh.H"
#include "vtkSurfaceWriter.H"
#include "meshTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const Foam::scalar Foam::internalAngleTolerance(80);
const Foam::scalar Foam::internalToleranceCosAngle
(
    cos(degToRad(180 - internalAngleTolerance))
);

const Foam::scalar Foam::externalAngleTolerance(10);
const Foam::scalar Foam::externalToleranceCosAngle
(
    cos(degToRad(180 - externalAngleTolerance))
);

void Foam::drawHitProblem
(
    const label fi,
    const triSurface& surf,
    const point& start,
    const point& p,
    const point& end,
    const pointIndexHitList& hitInfo
)
{
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

    meshTools::writeOBJ(Info, surf.points()[surf[fi][0]]);
    meshTools::writeOBJ(Info, surf.points()[surf[fi][1]]);
    meshTools::writeOBJ(Info, surf.points()[surf[fi][2]]);

    Info<< "f 4 5 6" << endl;

    forAll(hitInfo, hi)
    {
        label hFI = hitInfo[hi].index();

        meshTools::writeOBJ(Info, surf.points()[surf[hFI][0]]);
        meshTools::writeOBJ(Info, surf.points()[surf[hFI][1]]);
        meshTools::writeOBJ(Info, surf.points()[surf[hFI][2]]);

        Info<< "f "
            << 3*hi + 7 << " "
            << 3*hi + 8 << " "
            << 3*hi + 9
            << endl;
    }
}


void Foam::processHit
(
    scalar& internalCloseness,
    scalar& externalCloseness,
    const label fi,
    const triSurface& surf,
    const point& start,
    const point& p,
    const point& end,
    const vector& normal,
    const vectorField& normals,
    const pointIndexHitList& hitInfo
)
{
    if (hitInfo.size() < 1)
    {
        drawHitProblem(fi, surf, start, p, end, hitInfo);
    }
    else if (hitInfo.size() == 1)
    {
        if (!hitInfo[0].hit())
        {
        }
        else if (hitInfo[0].index() != fi)
        {
            drawHitProblem(fi, surf, start, p, end, hitInfo);
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
            drawHitProblem(fi, surf, start, p, end, hitInfo);
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


Foam::Pair<Foam::tmp<Foam::triSurfaceScalarField>> Foam::extractCloseness
(
    const triSurfaceMesh& surf
)
{
    const Time& runTime = surf.objectRegistry::time();

    // Prepare start and end points for intersection tests

    const vectorField& normals = surf.faceNormals();

    const scalar span = surf.bounds().mag();

    const pointField start(surf.faceCentres() - span*normals);
    const pointField end(surf.faceCentres() + span*normals);
    const pointField& faceCentres = surf.faceCentres();

    List<List<pointIndexHit>> allHitinfo;

    // Find all intersections (in order)
    surf.findLineAll(start, end, allHitinfo);

    scalarField internalCloseness(start.size(), great);
    scalarField externalCloseness(start.size(), great);

    forAll(allHitinfo, fi)
    {
        const List<pointIndexHit>& hitInfo = allHitinfo[fi];

        processHit
        (
            internalCloseness[fi],
            externalCloseness[fi],
            fi,
            surf,
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
                    surf.objectRegistry::name() + ".internalCloseness",
                    runTime.constant(),
                    "triSurface",
                    runTime
                ),
                surf,
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
                    surf.objectRegistry::name() + ".externalCloseness",
                    runTime.constant(),
                    "triSurface",
                    runTime
                ),
                surf,
                dimLength,
                externalCloseness
            )
        )
    );
}


// ************************************************************************* //
