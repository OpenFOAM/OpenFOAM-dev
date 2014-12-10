/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
    const char* NamedEnum<coupledPolyPatch::transformType, 5>::names[] =
    {
        "unknown",
        "rotational",
        "translational",
        "coincidentFullMatch",
        "noOrdering"
    };

    const NamedEnum<coupledPolyPatch::transformType, 5>
        coupledPolyPatch::transformTypeNames;
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

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
    const transformType transform
)
{
    pointField anchors(faces.size());

    if (transform != COINCIDENTFULLMATCH)
    {
        // Return the first point
        forAll(faces, faceI)
        {
            anchors[faceI] = points[faces[faceI][0]];
        }
    }
    else
    {
        // Make anchor point unique
        forAll(faces, faceI)
        {
            const face& f = faces[faceI];

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

                    // @todo Change to a tolerance and possibly select closest
                    // point to the origin
                    if (p1 == p2)
                    {
                        unique = false;
                        break;
                    }
                }

                if (unique)
                {
                    anchors[faceI] = p1;
                    break;
                }
            }

            if (!unique)
            {
                anchors[faceI] = points[faces[faceI][0]];
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

    forAll(faces, faceI)
    {
        const point& cc = faceCentres[faceI];

        const face& f = faces[faceI];

        // 1. calculate a typical size of the face. Use maximum distance
        //    to face centre
        scalar maxLenSqr = -GREAT;
        // 2. as measure of truncation error when comparing two coordinates
        //    use SMALL * maximum component
        scalar maxCmpt = -GREAT;

        forAll(f, fp)
        {
            const point& pt = points[f[fp]];
            maxLenSqr = max(maxLenSqr, magSqr(pt - cc));
            maxCmpt = max(maxCmpt, cmptMax(cmptMag(pt)));
        }

        tols[faceI] = max
        (
            SMALL,
            max(SMALL*maxCmpt, Foam::sqrt(maxLenSqr))
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
    scalar minDistSqr = GREAT;

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
                WarningIn
                (
                    "label coupledPolyPatch::getRotation\n"
                    "(\n"
                    "    const pointField&,\n"
                    "    const face&,\n"
                    "    const point&,\n"
                    "    const scalar\n"
                    ")"
                )   << "Cannot determine unique anchor point on face "
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


void Foam::coupledPolyPatch::calcTransformTensors
(
    const vectorField& Cf,
    const vectorField& Cr,
    const vectorField& nf,
    const vectorField& nr,
    const scalarField& smallDist,
    const scalar absTol,
    const transformType transform
) const
{
    if (debug)
    {
        Pout<< "coupledPolyPatch::calcTransformTensors : " << name() << endl
            << "    transform:" << transformTypeNames[transform] << nl
            << "    (half)size:" << Cf.size() << nl
            << "    absTol:" << absTol << nl
            << "    smallDist min:" << min(smallDist) << nl
            << "    smallDist max:" << max(smallDist) << nl
            << "    sum(mag(nf & nr)):" << sum(mag(nf & nr)) << endl;
    }

    // Tolerance calculation.
    // - normal calculation: assume absTol is the absolute error in a
    // single normal/transformation calculation. Consists both of numerical
    // precision (on the order of SMALL and of writing precision
    // (from e.g. decomposition)
    // Then the overall error of summing the normals is sqrt(size())*absTol
    // - separation calculation: pass in from the outside an allowable error.

    if (Cf.size() == 0)
    {
        // Dummy geometry. Assume non-separated, parallel.
        separation_.setSize(0);
        forwardT_.clear();
        reverseT_.clear();
        collocated_.setSize(0);
    }
    else
    {
        scalar error = absTol*Foam::sqrt(1.0*Cf.size());

        if (debug)
        {
            Pout<< "    error:" << error << endl;
        }

        if
        (
            transform == ROTATIONAL
         || (
                transform != TRANSLATIONAL
             && transform != COINCIDENTFULLMATCH
             && (sum(mag(nf & nr)) < Cf.size() - error)
            )
        )
        {
            // Type is rotation or unknown and normals not aligned

            // Assume per-face differing transformation, correct later

            separation_.setSize(0);

            forwardT_.setSize(Cf.size());
            reverseT_.setSize(Cf.size());
            collocated_.setSize(Cf.size());
            collocated_ = false;

            forAll(forwardT_, facei)
            {
                forwardT_[facei] = rotationTensor(-nr[facei], nf[facei]);
                reverseT_[facei] = rotationTensor(nf[facei], -nr[facei]);
            }

            if (debug)
            {
                Pout<< "    sum(mag(forwardT_ - forwardT_[0])):"
                    << sum(mag(forwardT_ - forwardT_[0]))
                    << endl;
            }

            if (sum(mag(forwardT_ - forwardT_[0])) < error)
            {
                forwardT_.setSize(1);
                reverseT_.setSize(1);
                collocated_.setSize(1);

                if (debug)
                {
                    Pout<< "    difference in rotation less than"
                        << " local tolerance "
                        << error << ". Assuming uniform rotation." << endl;
                }
            }
        }
        else
        {
            // Translational or (unknown and normals aligned)

            forwardT_.setSize(0);
            reverseT_.setSize(0);

            separation_ = Cr - Cf;

            collocated_.setSize(separation_.size());

            // Three situations:
            // - separation is zero. No separation.
            // - separation is same. Single separation vector.
            // - separation differs per face. Separation vectorField.

            // Check for different separation per face
            bool sameSeparation = true;
            bool doneWarning = false;

            forAll(separation_, facei)
            {
                scalar smallSqr = sqr(smallDist[facei]);

                collocated_[facei] = (magSqr(separation_[facei]) < smallSqr);

                // Check if separation differing w.r.t. face 0.
                if (magSqr(separation_[facei] - separation_[0]) > smallSqr)
                {
                    sameSeparation = false;

                    if (!doneWarning && debug)
                    {
                        doneWarning = true;

                        Pout<< "    separation " << separation_[facei]
                            << " at " << facei
                            << " differs from separation[0] " << separation_[0]
                            << " by more than local tolerance "
                            << smallDist[facei]
                            << ". Assuming non-uniform separation." << endl;
                    }
                }
            }

            if (sameSeparation)
            {
                // Check for zero separation (at 0 so everywhere)
                if (collocated_[0])
                {
                    if (debug)
                    {
                        Pout<< "    separation " << mag(separation_[0])
                            << " less than local tolerance " << smallDist[0]
                            << ". Assuming zero separation." << endl;
                    }

                    separation_.setSize(0);
                    collocated_ = boolList(1, true);
                }
                else
                {
                    if (debug)
                    {
                        Pout<< "    separation " << mag(separation_[0])
                            << " more than local tolerance " << smallDist[0]
                            << ". Assuming uniform separation." << endl;
                    }

                    separation_.setSize(1);
                    collocated_ = boolList(1, false);
                }
            }
        }
    }

    if (debug)
    {
        Pout<< "    separation_:" << separation_.size() << nl
            << "    forwardT size:" << forwardT_.size() << endl;
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
    const transformType transform
)
:
    polyPatch(name, size, start, index, bm, patchType),
    matchTolerance_(defaultMatchTol_),
    transform_(transform)
{}


Foam::coupledPolyPatch::coupledPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    polyPatch(name, dict, index, bm, patchType),
    matchTolerance_(dict.lookupOrDefault("matchTolerance", defaultMatchTol_)),
    transform_
    (
        dict.found("transform")
      ? transformTypeNames.read(dict.lookup("transform"))
      : UNKNOWN
    )
{}


Foam::coupledPolyPatch::coupledPolyPatch
(
    const coupledPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    polyPatch(pp, bm),
    matchTolerance_(pp.matchTolerance_),
    transform_(pp.transform_)
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
    matchTolerance_(pp.matchTolerance_),
    transform_(pp.transform_)
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
    matchTolerance_(pp.matchTolerance_),
    transform_(pp.transform_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coupledPolyPatch::~coupledPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::coupledPolyPatch::write(Ostream& os) const
{
    polyPatch::write(os);
    //if (matchTolerance_ != defaultMatchTol_)
    {
        os.writeKeyword("matchTolerance") << matchTolerance_
            << token::END_STATEMENT << nl;
        os.writeKeyword("transform") << transformTypeNames[transform_]
            << token::END_STATEMENT << nl;
    }
}


// ************************************************************************* //
