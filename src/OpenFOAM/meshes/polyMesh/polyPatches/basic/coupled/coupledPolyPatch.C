/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "coupledPolyPatch.H"
#include "ListOps.H"
#include "transform.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(coupledPolyPatch, 0);

    const scalar coupledPolyPatch::defaultMatchTol_ = 1e-4;

    template<>
    const char* NamedEnum<coupledPolyPatch::orderingType, 3>::names[] =
    {
        "unknown",
        "coincidentFullMatch",
        "noOrdering"
    };

    const NamedEnum<coupledPolyPatch::orderingType, 3>
        coupledPolyPatch::orderingTypeNames;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::coupledPolyPatch::writeOBJ(Ostream& os, const point& pt)
{
    os << "v " << pt.x() << ' ' << pt.y() << ' ' << pt.z() << endl;
}


void Foam::coupledPolyPatch::writeOBJ
(
    Ostream& os,
    const pointField& points,
    const labelList& pointLabels
)
{
    forAll(pointLabels, i)
    {
        writeOBJ(os, points[pointLabels[i]]);
    }
}


void Foam::coupledPolyPatch::writeOBJ
(
    Ostream& os,
    const point& p0,
    const point& p1,
    label& vertI
)
{
    writeOBJ(os, p0);
    vertI++;

    writeOBJ(os, p1);
    vertI++;

    os  << "l " << vertI-1 << ' ' << vertI << nl;
}


void Foam::coupledPolyPatch::writeOBJ
(
    const fileName& fName,
    const UList<face>& faces,
    const pointField& points
)
{
    OFstream os(fName);

    Map<label> foamToObj(4*faces.size());

    label vertI = 0;

    forAll(faces, i)
    {
        const face& f = faces[i];

        forAll(f, fp)
        {
            if (foamToObj.insert(f[fp], vertI))
            {
                writeOBJ(os, points[f[fp]]);
                vertI++;
            }
        }

        os << 'l';
        forAll(f, fp)
        {
            os << ' ' << foamToObj[f[fp]]+1;
        }
        os << ' ' << foamToObj[f[0]]+1 << nl;
    }
}


Foam::pointField Foam::coupledPolyPatch::getAnchorPoints
(
    const UList<face>& faces,
    const pointField& points,
    const orderingType ordering
)
{
    pointField anchors(faces.size());

    if (ordering != COINCIDENTFULLMATCH)
    {
        // Return the first point
        forAll(faces, facei)
        {
            anchors[facei] = points[faces[facei][0]];
        }
    }
    else
    {
        // Make anchor point unique
        forAll(faces, facei)
        {
            const face& f = faces[facei];

            bool unique = true;

            forAll(f, fp1)
            {
                const point& p1 = points[f[fp1]];

                unique = true;

                for (label fp2 = 0; fp2 < f.size(); ++fp2)
                {
                    if (f[fp1] == f[fp2])
                    {
                        continue;
                    }

                    const point& p2 = points[f[fp2]];

                    // TODO: Change to a tolerance and possibly select closest
                    // point to the origin
                    if (p1 == p2)
                    {
                        unique = false;
                        break;
                    }
                }

                if (unique)
                {
                    anchors[facei] = p1;
                    break;
                }
            }

            if (!unique)
            {
                anchors[facei] = points[faces[facei][0]];
            }
        }
    }

    return anchors;
}


Foam::scalarField Foam::coupledPolyPatch::calcFaceTol
(
    const UList<face>& faces,
    const pointField& points,
    const pointField& faceCentres
)
{
    // Calculate typical distance per face
    scalarField tols(faces.size());

    forAll(faces, facei)
    {
        const point& cc = faceCentres[facei];

        const face& f = faces[facei];

        // 1. calculate a typical size of the face. Use maximum distance
        //    to face centre
        scalar maxLenSqr = -great;
        // 2. as measure of truncation error when comparing two coordinates
        //    use small * maximum component
        scalar maxCmpt = -great;

        forAll(f, fp)
        {
            const point& pt = points[f[fp]];
            maxLenSqr = max(maxLenSqr, magSqr(pt - cc));
            maxCmpt = max(maxCmpt, cmptMax(cmptMag(pt)));
        }

        tols[facei] = max
        (
            small,
            max(small*maxCmpt, Foam::sqrt(maxLenSqr))
        );
    }
    return tols;
}


Foam::label Foam::coupledPolyPatch::getRotation
(
    const pointField& points,
    const face& f,
    const point& anchor,
    const scalar tol
)
{
    label anchorFp = -1;
    scalar minDistSqr = great;

    forAll(f, fp)
    {
        scalar distSqr = magSqr(anchor - points[f[fp]]);

        if (distSqr < minDistSqr)
        {
            minDistSqr = distSqr;
            anchorFp = fp;
        }
    }

    if (anchorFp == -1 || Foam::sqrt(minDistSqr) > tol)
    {
        return -1;
    }
    else
    {
        // Check that anchor is unique.
        forAll(f, fp)
        {
            scalar distSqr = magSqr(anchor - points[f[fp]]);

            if (distSqr == minDistSqr && fp != anchorFp)
            {
                WarningInFunction
                    << "Cannot determine unique anchor point on face "
                    << UIndirectList<point>(points, f)
                    << endl
                    << "Both at index " << anchorFp << " and " << fp
                    << " the vertices have the same distance "
                    << Foam::sqrt(minDistSqr)
                    << " to the anchor " << anchor
                    << ". Continuing but results might be wrong."
                    << nl << endl;
            }
        }

        // Positive rotation
        return (f.size() - anchorFp) % f.size();
    }
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::coupledPolyPatch::coupledPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType,
    const orderingType ordering
)
:
    polyPatch(name, size, start, index, bm, patchType),
    matchTolerance_(defaultMatchTol_),
    ordering_(ordering)
{}


Foam::coupledPolyPatch::coupledPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType,
    const orderingType ordering
)
:
    polyPatch(name, dict, index, bm, patchType),
    matchTolerance_(dict.lookupOrDefault("matchTolerance", defaultMatchTol_)),
    ordering_
    (
        dict.found("ordering")
      ? orderingTypeNames.read(dict.lookup("ordering"))
      : ordering
    )
{}


Foam::coupledPolyPatch::coupledPolyPatch
(
    const coupledPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    polyPatch(pp, bm),
    matchTolerance_(pp.matchTolerance_)
{}


Foam::coupledPolyPatch::coupledPolyPatch
(
    const coupledPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    polyPatch(pp, bm, index, newSize, newStart),
    matchTolerance_(pp.matchTolerance_)
{}


Foam::coupledPolyPatch::coupledPolyPatch
(
    const coupledPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    polyPatch(pp, bm, index, mapAddressing, newStart),
    matchTolerance_(pp.matchTolerance_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coupledPolyPatch::~coupledPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::coupledPolyPatch::write(Ostream& os) const
{
    polyPatch::write(os);
    writeEntry(os, "matchTolerance", matchTolerance_);
    writeEntry(os, "ordering", orderingTypeNames[ordering_]);
}


// ************************************************************************* //
