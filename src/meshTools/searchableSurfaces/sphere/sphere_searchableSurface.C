/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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

#include "sphere_searchableSurface.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace searchableSurfaces
    {
        defineTypeNameAndDebug(sphere, 0);

        addToRunTimeSelectionTable
        (
            searchableSurface,
            sphere,
            dictionary
        );

        addBackwardCompatibleToRunTimeSelectionTable
        (
            searchableSurface,
            sphere,
            dictionary,
            searchableSphere,
            "searchableSphere"
        );
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::pointIndexHit Foam::searchableSurfaces::sphere::findNearest
(
    const point& sample,
    const scalar nearestDistSqr
) const
{
    pointIndexHit info(false, sample, -1);

    const vector n(sample - centre_);
    scalar magN = mag(n);

    if (nearestDistSqr >= sqr(magN - radius_))
    {
        if (magN < rootVSmall)
        {
            info.rawPoint() = centre_ + vector(1,0,0)*radius_;
        }
        else
        {
            info.rawPoint() = centre_ + n/magN*radius_;
        }
        info.setHit();
        info.setIndex(0);
    }

    return info;
}


// From Graphics Gems - intersection of sphere with ray
void Foam::searchableSurfaces::sphere::findLineAll
(
    const point& start,
    const point& end,
    pointIndexHit& near,
    pointIndexHit& far
) const
{
    near.setMiss();
    far.setMiss();

    vector dir(end-start);
    scalar magSqrDir = magSqr(dir);

    if (magSqrDir > rootVSmall)
    {
        const vector toCentre(centre_-start);
        scalar magSqrToCentre = magSqr(toCentre);

        dir /= Foam::sqrt(magSqrDir);

        scalar v = (toCentre & dir);

        scalar disc = sqr(radius_) - (magSqrToCentre - sqr(v));

        if (disc >= 0)
        {
            scalar d = Foam::sqrt(disc);

            scalar nearParam = v-d;

            if (nearParam >= 0 && sqr(nearParam) <= magSqrDir)
            {
                near.setHit();
                near.setPoint(start + nearParam*dir);
                near.setIndex(0);
            }

            scalar farParam = v+d;

            if (farParam >= 0 && sqr(farParam) <= magSqrDir)
            {
                far.setHit();
                far.setPoint(start + farParam*dir);
                far.setIndex(0);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::searchableSurfaces::sphere::sphere
(
    const IOobject& io,
    const point& centre,
    const scalar radius
)
:
    searchableSurface(io),
    centre_(centre),
    radius_(radius)
{
    bounds() = boundBox
    (
        centre_ - radius_*vector::one,
        centre_ + radius_*vector::one
    );
}


Foam::searchableSurfaces::sphere::sphere
(
    const IOobject& io,
    const dictionary& dict
)
:
    searchableSurface(io),
    centre_(dict.lookup<point>("centre", dimLength)),
    radius_(dict.lookup<scalar>("radius", dimLength))
{
    bounds() = boundBox
    (
        centre_ - radius_*vector::one,
        centre_ + radius_*vector::one
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::searchableSurfaces::sphere::~sphere()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::searchableSurfaces::sphere::overlaps(const boundBox& bb) const
{
    return bb.overlaps(centre_, sqr(radius_));
}


const Foam::wordList& Foam::searchableSurfaces::sphere::regions() const
{
    if (regions_.empty())
    {
        regions_.setSize(1);
        regions_[0] = "region0";
    }
    return regions_;
}



void Foam::searchableSurfaces::sphere::boundingSpheres
(
    pointField& centres,
    scalarField& radiusSqr
) const
{
    centres.setSize(1);
    centres[0] = centre_;

    radiusSqr.setSize(1);
    radiusSqr[0] = Foam::sqr(radius_);

    // Add a bit to make sure all points are tested inside
    radiusSqr += Foam::sqr(small);
}


void Foam::searchableSurfaces::sphere::findNearest
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    List<pointIndexHit>& info
) const
{
    info.setSize(samples.size());

    forAll(samples, i)
    {
        info[i] = findNearest(samples[i], nearestDistSqr[i]);
    }
}


void Foam::searchableSurfaces::sphere::findLine
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    info.setSize(start.size());

    pointIndexHit b;

    forAll(start, i)
    {
        // Pick nearest intersection. If none intersected take second one.
        findLineAll(start[i], end[i], info[i], b);
        if (!info[i].hit() && b.hit())
        {
            info[i] = b;
        }
    }
}


void Foam::searchableSurfaces::sphere::findLineAny
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    info.setSize(start.size());

    pointIndexHit b;

    forAll(start, i)
    {
        // Discard far intersection
        findLineAll(start[i], end[i], info[i], b);
        if (!info[i].hit() && b.hit())
        {
            info[i] = b;
        }
    }
}


void Foam::searchableSurfaces::sphere::findLineAll
(
    const pointField& start,
    const pointField& end,
    List<List<pointIndexHit>>& info
) const
{
    info.setSize(start.size());

    forAll(start, i)
    {
        pointIndexHit near, far;
        findLineAll(start[i], end[i], near, far);

        if (near.hit())
        {
            if (far.hit())
            {
                info[i].setSize(2);
                info[i][0] = near;
                info[i][1] = far;
            }
            else
            {
                info[i].setSize(1);
                info[i][0] = near;
            }
        }
        else
        {
            if (far.hit())
            {
                info[i].setSize(1);
                info[i][0] = far;
            }
            else
            {
                info[i].clear();
            }
        }
    }
}


void Foam::searchableSurfaces::sphere::getRegion
(
    const List<pointIndexHit>& info,
    labelList& region
) const
{
    region.setSize(info.size());
    region = 0;
}


void Foam::searchableSurfaces::sphere::getNormal
(
    const List<pointIndexHit>& info,
    vectorField& normal
) const
{
    normal.setSize(info.size());
    normal = Zero;

    forAll(info, i)
    {
        if (info[i].hit())
        {
            normal[i] = info[i].hitPoint() - centre_;

            normal[i] /= mag(normal[i])+vSmall;
        }
        else
        {
            // Set to what?
        }
    }
}


void Foam::searchableSurfaces::sphere::getVolumeType
(
    const pointField& points,
    List<volumeType>& volType
) const
{
    volType.setSize(points.size());
    volType = volumeType::inside;

    forAll(points, pointi)
    {
        const point& pt = points[pointi];

        if (magSqr(pt - centre_) <= sqr(radius_))
        {
            volType[pointi] = volumeType::inside;
        }
        else
        {
            volType[pointi] = volumeType::outside;
        }
    }
}


// ************************************************************************* //
