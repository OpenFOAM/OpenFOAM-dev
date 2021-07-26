/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "surfaceToCell.H"
#include "polyMesh.H"
#include "meshSearch.H"
#include "triSurface.H"
#include "triSurfaceSearch.H"
#include "cellClassification.H"
#include "cpuTime.H"
#include "demandDrivenData.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(surfaceToCell, 0);
    addToRunTimeSelectionTable(topoSetSource, surfaceToCell, word);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::surfaceToCell::getNearest
(
    const triSurfaceSearch& querySurf,
    const label pointi,
    const point& pt,
    const vector& span,
    Map<label>& cache
)
{
    Map<label>::const_iterator iter = cache.find(pointi);

    if (iter != cache.end())
    {
        // Found cached answer
        return iter();
    }
    else
    {
        pointIndexHit inter = querySurf.nearest(pt, span);

        // Triangle label (can be -1)
        label triI = inter.index();

        // Store triangle on point
        cache.insert(pointi, triI);

        return triI;
    }
}


bool Foam::surfaceToCell::differingPointNormals
(
    const triSurfaceSearch& querySurf,

    const vector& span,         // Current search span
    const label celli,
    const label cellTriI,       // Nearest (to cell centre) surface triangle

    Map<label>& pointToNearest  // Cache for nearest triangle to point
) const
{
    const triSurface& surf = querySurf.surface();
    const vectorField& normals = surf.faceNormals();

    const faceList& faces = mesh().faces();
    const pointField& points = mesh().points();

    const labelList& cFaces = mesh().cells()[celli];

    forAll(cFaces, cFacei)
    {
        const face& f = faces[cFaces[cFacei]];

        forAll(f, fp)
        {
            label pointi = f[fp];

            label pointTriI =
                getNearest
                (
                    querySurf,
                    pointi,
                    points[pointi],
                    span,
                    pointToNearest
                );

            if (pointTriI != -1 && pointTriI != cellTriI)
            {
                scalar cosAngle = normals[pointTriI] & normals[cellTriI];

                if (cosAngle < 0.9)
                {
                    return true;
                }
            }
        }
    }
    return false;
}


void Foam::surfaceToCell::combine(topoSet& set, const bool add) const
{
    cpuTime timer;


    if (useSurfaceOrientation_ && (includeInside_ || includeOutside_))
    {
        const meshSearch queryMesh(mesh_);

        //- Calculate for each searchPoint inside/outside status.
        boolList isInside(querySurf().calcInside(mesh_.cellCentres()));

        Info<< "    Marked inside/outside using surface orientation in = "
            << timer.cpuTimeIncrement() << " s" << endl << endl;

        forAll(isInside, celli)
        {
            if (isInside[celli] && includeInside_)
            {
                addOrDelete(set, celli, add);
            }
            else if (!isInside[celli] && includeOutside_)
            {
                addOrDelete(set, celli, add);
            }
        }
    }
    else if (includeCut_ || includeInside_ || includeOutside_)
    {
        //
        // Cut cells with surface and classify cells
        //


        // Construct search engine on mesh

        const meshSearch queryMesh(mesh_);


        // Check all 'outside' points
        forAll(outsidePoints_, outsideI)
        {
            const point& outsidePoint = outsidePoints_[outsideI];

            // Find cell point is in. Linear search.
            label celli = queryMesh.findCell(outsidePoint, -1, false);
            if (returnReduce(celli, maxOp<label>()) == -1)
            {
                FatalErrorInFunction
                    << "outsidePoint " << outsidePoint
                    << " is not inside any cell"
                    << exit(FatalError);
            }
        }

        // Cut faces with surface and classify cells

        cellClassification cellType
        (
            mesh_,
            queryMesh,
            querySurf(),
            outsidePoints_
        );


        Info<< "    Marked inside/outside using surface intersection in = "
            << timer.cpuTimeIncrement() << " s" << endl << endl;

        //- Add/remove cells using set
        forAll(cellType, celli)
        {
            label cType = cellType[celli];

            if
            (
                (
                    includeCut_
                 && (cType == cellClassification::CUT)
                )
             || (
                    includeInside_
                 && (cType == cellClassification::INSIDE)
                )
             || (
                    includeOutside_
                 && (cType == cellClassification::OUTSIDE)
                )
            )
            {
                addOrDelete(set, celli, add);
            }
        }
    }


    if (nearDist_ > 0)
    {
        //
        // Determine distance to surface
        //

        const pointField& ctrs = mesh_.cellCentres();

        // Box dimensions to search in octree.
        const vector span(nearDist_, nearDist_, nearDist_);


        if (curvature_ < -1)
        {
            Info<< "    Selecting cells with cellCentre closer than "
                << nearDist_ << " to surface" << endl;

            // No need to test curvature. Insert near cells into set.

            forAll(ctrs, celli)
            {
                const point& c = ctrs[celli];

                pointIndexHit inter = querySurf().nearest(c, span);

                if (inter.hit() && (mag(inter.hitPoint() - c) < nearDist_))
                {
                    addOrDelete(set, celli, add);
                }
            }

            Info<< "    Determined nearest surface point in = "
                << timer.cpuTimeIncrement() << " s" << endl << endl;

        }
        else
        {
            // Test near cells for curvature

            Info<< "    Selecting cells with cellCentre closer than "
                << nearDist_ << " to surface and curvature factor"
                << " less than " << curvature_ << endl;

            // Cache for nearest surface triangle for a point
            Map<label> pointToNearest(mesh_.nCells()/10);

            forAll(ctrs, celli)
            {
                const point& c = ctrs[celli];

                pointIndexHit inter = querySurf().nearest(c, span);

                if (inter.hit() && (mag(inter.hitPoint() - c) < nearDist_))
                {
                    if
                    (
                        differingPointNormals
                        (
                            querySurf(),
                            span,
                            celli,
                            inter.index(),      // nearest surface triangle
                            pointToNearest
                        )
                    )
                    {
                        addOrDelete(set, celli, add);
                    }
                }
            }

            Info<< "    Determined nearest surface point in = "
                << timer.cpuTimeIncrement() << " s" << endl << endl;
        }
    }
}


void Foam::surfaceToCell::checkSettings() const
{
    if
    (
        (nearDist_ < 0)
     && (curvature_ < -1)
     && (
            (includeCut_ && includeInside_ && includeOutside_)
         || (!includeCut_ && !includeInside_ && !includeOutside_)
        )
    )
    {
        FatalErrorInFunction
            << "Illegal include cell specification."
            << " Result would be either all or no cells." << endl
            << "Please set one of includeCut, includeInside, includeOutside"
            << " to true, set nearDistance to a value > 0"
            << " or set curvature to a value -1 .. 1."
            << exit(FatalError);
    }

    if (useSurfaceOrientation_ && includeCut_)
    {
        FatalErrorInFunction
            << "Illegal include cell specification."
            << " You cannot specify both 'useSurfaceOrientation'"
            << " and 'includeCut'"
            << " since 'includeCut' specifies a topological split"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceToCell::surfaceToCell
(
    const polyMesh& mesh,
    const fileName& surfName,
    const pointField& outsidePoints,
    const bool includeCut,
    const bool includeInside,
    const bool includeOutside,
    const bool useSurfaceOrientation,
    const scalar nearDist,
    const scalar curvature
)
:
    topoSetSource(mesh),
    surfName_(surfName),
    outsidePoints_(outsidePoints),
    includeCut_(includeCut),
    includeInside_(includeInside),
    includeOutside_(includeOutside),
    useSurfaceOrientation_(useSurfaceOrientation),
    nearDist_(nearDist),
    curvature_(curvature),
    surfPtr_(new triSurface(surfName_)),
    querySurfPtr_(new triSurfaceSearch(*surfPtr_)),
    IOwnPtrs_(true)
{
    checkSettings();
}


Foam::surfaceToCell::surfaceToCell
(
    const polyMesh& mesh,
    const fileName& surfName,
    const triSurface& surf,
    const triSurfaceSearch& querySurf,
    const pointField& outsidePoints,
    const bool includeCut,
    const bool includeInside,
    const bool includeOutside,
    const bool useSurfaceOrientation,
    const scalar nearDist,
    const scalar curvature
)
:
    topoSetSource(mesh),
    surfName_(surfName),
    outsidePoints_(outsidePoints),
    includeCut_(includeCut),
    includeInside_(includeInside),
    includeOutside_(includeOutside),
    useSurfaceOrientation_(useSurfaceOrientation),
    nearDist_(nearDist),
    curvature_(curvature),
    surfPtr_(&surf),
    querySurfPtr_(&querySurf),
    IOwnPtrs_(false)
{
    checkSettings();
}


Foam::surfaceToCell::surfaceToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    surfName_(fileName(dict.lookup("file")).expand()),
    outsidePoints_(dict.lookup("outsidePoints")),
    includeCut_(readBool(dict.lookup("includeCut"))),
    includeInside_(readBool(dict.lookup("includeInside"))),
    includeOutside_(readBool(dict.lookup("includeOutside"))),
    useSurfaceOrientation_
    (
        dict.lookupOrDefault<bool>("useSurfaceOrientation", false)
    ),
    nearDist_(dict.lookup<scalar>("nearDistance")),
    curvature_(dict.lookup<scalar>("curvature")),
    surfPtr_(new triSurface(surfName_)),
    querySurfPtr_(new triSurfaceSearch(*surfPtr_)),
    IOwnPtrs_(true)
{
    checkSettings();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfaceToCell::~surfaceToCell()
{
    if (IOwnPtrs_)
    {
        deleteDemandDrivenData(surfPtr_);
        deleteDemandDrivenData(querySurfPtr_);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::surfaceToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding cells in relation to surface " << surfName_
            << " ..." << endl;

        combine(set, true);
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing cells in relation to surface " << surfName_
            << " ..." << endl;

        combine(set, false);
    }
}


// ************************************************************************* //
