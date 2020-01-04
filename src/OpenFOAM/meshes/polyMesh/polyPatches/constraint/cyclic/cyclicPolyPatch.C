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

#include "cyclicPolyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "polyBoundaryMesh.H"
#include "polyMesh.H"
#include "demandDrivenData.H"
#include "OFstream.H"
#include "matchPoints.H"
#include "EdgeMap.H"
#include "Time.H"
#include "transformField.H"
#include "SubField.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, cyclicPolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, cyclicPolyPatch, dictionary);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::cyclicPolyPatch::calcTransformTensors
(
    const vectorField& thisPatchCtrs,
    const vectorField& nbrPatchCtrs,
    const vectorField& thisPatchNormals,
    const vectorField& nbrPatchNormals,
    const scalarField& smallDist,
    const scalar absTol,
    const orderingType ordering,
    const transformTypes transform
) const
{
    if (debug)
    {
        Pout<< "coupledPolyPatch::calcTransformTensors : " << name() << endl
            << "    transform:" << transformTypeNames[transform] << nl
            << "    (half)size:" << thisPatchCtrs.size() << nl
            << "    absTol:" << absTol << nl
            << "    smallDist min:" << min(smallDist) << nl
            << "    smallDist max:" << max(smallDist) << nl
            << "    sum(mag(thisPatchNormals & nbrPatchNormals)):"
            << sum(mag(thisPatchNormals & nbrPatchNormals)) << endl;
    }

    // Tolerance calculation.
    // - normal calculation: assume absTol is the absolute error in a
    // single normal/transformation calculation. Consists both of numerical
    // precision (on the order of small and of writing precision
    // (from e.g. decomposition)
    // Then the overall error of summing the normals is sqrt(size())*absTol
    // - separation calculation: pass in from the outside an allowable error.

    if (thisPatchCtrs.size() == 0)
    {
        // Dummy geometry. Assume non-separated, parallel.
        transform_ = transformer();
    }
    else
    {
        scalar error = absTol*Foam::sqrt(1.0*thisPatchCtrs.size());

        if (debug)
        {
            Pout<< "    error:" << error << endl;
        }

        if
        (
            transform == ROTATIONAL
         || (
                transform != TRANSLATIONAL
             && ordering != COINCIDENTFULLMATCH
             && (
                    sum(mag(thisPatchNormals & nbrPatchNormals))
                  < thisPatchCtrs.size() - error
                )
            )
        )
        {
            // Type is rotation or unknown and normals not aligned

            tensorField forwardT(thisPatchCtrs.size());
            tensorField reverseT(thisPatchCtrs.size());

            forAll(forwardT, facei)
            {
                forwardT[facei] = rotationTensor
                (
                    -nbrPatchNormals[facei],
                    thisPatchNormals[facei]
                );

                reverseT[facei] = rotationTensor
                (
                    thisPatchNormals[facei],
                    -nbrPatchNormals[facei]
                );
            }

            if (sum(mag(forwardT - forwardT[0])) > error)
            {
                Pout<< "--> FOAM Warning : "
                    << " Variation in rotation greater than"
                    << " local tolerance " << error << endl;
            }

            transform_ = transformer(forwardT[0]);
        }
        else
        {
            // Translational or (unknown and normals aligned)

            // Three situations:
            // - separation is zero. No separation.
            // - separation is same. Single separation vector.
            // - separation differs per face -> error.

            // Check for different separation per face
            bool sameSeparation = true;
            bool doneWarning = false;

            const vectorField separation(nbrPatchCtrs - thisPatchCtrs);

            forAll(separation, facei)
            {
                const scalar smallSqr = sqr(smallDist[facei]);

                // Check if separation differing w.r.t. face 0.
                if (magSqr(separation[facei] - separation[0]) > smallSqr)
                {
                    sameSeparation = false;

                    if (!doneWarning && debug)
                    {
                        doneWarning = true;

                        Pout<< "    separation " << separation[facei]
                            << " at " << facei
                            << " differs from separation[0] " << separation[0]
                            << " by more than local tolerance "
                            << smallDist[facei]
                            << ". Assuming non-uniform separation." << endl;
                    }
                }
            }

            if (sameSeparation)
            {
                // Check for zero separation
                if (mag(separation[0]) < smallDist[0])
                {
                    if (debug)
                    {
                        Pout<< "    separation " << mag(separation[0])
                            << " less than local tolerance " << smallDist[0]
                            << ". Assuming zero separation." << endl;
                    }

                    transform_ = transformer();
                }
                else
                {
                    if (debug)
                    {
                        Pout<< "    separation " << mag(separation[0])
                            << " more than local tolerance " << smallDist[0]
                            << ". Assuming uniform separation." << endl;
                    }

                    transform_ = transformer(separation[0]);
                }
            }
            else
            {
                Pout<< "--> FOAM Warning : "
                    << " Variation in separation greater than"
                    << " local tolerance " << smallDist[0] << endl;

                transform_ = transformer(separation[0]);
            }
        }
    }
}


Foam::label Foam::cyclicPolyPatch::findMaxArea
(
    const pointField& points,
    const faceList& faces
)
{
    label maxI = -1;
    scalar maxAreaSqr = -great;

    forAll(faces, facei)
    {
        scalar areaSqr = magSqr(faces[facei].area(points));

        if (areaSqr > maxAreaSqr)
        {
            maxAreaSqr = areaSqr;
            maxI = facei;
        }
    }
    return maxI;
}


void Foam::cyclicPolyPatch::calcTransforms()
{
    if (size())
    {
        // thisPatch
        const cyclicPolyPatch& thisPatch = *this;
        vectorField thisPatchAreas(thisPatch.size());
        forAll(thisPatch, facei)
        {
            thisPatchAreas[facei] = thisPatch[facei].area(thisPatch.points());
        }

        // nbrPatch
        const cyclicPolyPatch& nbrPatch = this->nbrPatch();
        vectorField nbrPatchAreas(nbrPatch.size());
        forAll(nbrPatch, facei)
        {
            nbrPatchAreas[facei] = nbrPatch[facei].area(nbrPatch.points());
        }

        calcTransforms
        (
            thisPatch,
            thisPatch.faceCentres(),
            thisPatchAreas,
            nbrPatch.faceCentres(),
            nbrPatchAreas
        );
    }
}


void Foam::cyclicPolyPatch::calcTransforms
(
    const primitivePatch& thisPatch,
    const pointField& thisPatchCtrs,
    const vectorField& thisPatchAreas,
    const pointField& nbrPatchCtrs,
    const vectorField& nbrPatchAreas
)
{
    if (debug && owner())
    {
        fileName casePath(boundaryMesh().mesh().time().path());
        {
            fileName nm0(casePath/name()+"_faces.obj");
            Pout<< "cyclicPolyPatch::calcTransforms : Writing " << name()
                << " faces to OBJ file " << nm0 << endl;
            writeOBJ(nm0, thisPatch, thisPatch.points());
        }
        const cyclicPolyPatch& nbrPatch = this->nbrPatch();
        {
            fileName nm1(casePath/nbrPatch.name()+"_faces.obj");
            Pout<< "cyclicPolyPatch::calcTransforms : Writing "
                << nbrPatch.name()
                << " faces to OBJ file " << nm1 << endl;
            writeOBJ(nm1, nbrPatch, nbrPatch.points());
        }
        {
            OFstream str(casePath/name()+"_to_" + nbrPatch.name() + ".obj");
            label vertI = 0;
            Pout<< "cyclicPolyPatch::calcTransforms :"
                << " Writing coupled face centres as lines to " << str.name()
                << endl;
            forAll(thisPatchCtrs, i)
            {
                const point& p0 = thisPatchCtrs[i];
                str << "v " << p0.x() << ' ' << p0.y() << ' ' << p0.z() << nl;
                vertI++;
                const point& p1 = nbrPatchCtrs[i];
                str << "v " << p1.x() << ' ' << p1.y() << ' ' << p1.z() << nl;
                vertI++;
                str << "l " << vertI-1 << ' ' << vertI << nl;
            }
        }
    }


    // Some sanity checks

    if (thisPatchCtrs.size() != nbrPatchCtrs.size())
    {
        FatalErrorInFunction
            << "For patch " << name()
            << " there are " << thisPatchCtrs.size()
            << " face centres, for the neighbour patch " << nbrPatch().name()
            << " there are " << nbrPatchCtrs.size()
            << exit(FatalError);
    }

    if (transformType() != nbrPatch().transformType())
    {
        FatalErrorInFunction
            << "Patch " << name()
            << " has transform type " << transformTypeNames[transformType()]
            << ", neighbour patch " << nbrPatchName()
            << " has transform type "
            << nbrPatch().transformTypeNames[nbrPatch().transformType()]
            << exit(FatalError);
    }


    // Calculate transformation tensors

    if (thisPatchCtrs.size() > 0)
    {
        vectorField thisPatchNormals(thisPatchAreas.size());
        vectorField nbrPatchNormals(nbrPatchAreas.size());

        scalar maxAreaDiff = -great;
        label maxAreaFacei = -1;

        forAll(thisPatch, facei)
        {
            scalar magSf = mag(thisPatchAreas[facei]);
            scalar nbrMagSf = mag(nbrPatchAreas[facei]);
            scalar avSf = (magSf + nbrMagSf)/2.0;

            if (magSf < rootVSmall && nbrMagSf < rootVSmall)
            {
                // Undetermined normal. Use dummy normal to force separation
                // check. (note use of sqrt(vSmall) since that is how mag
                // scales)
                thisPatchNormals[facei] = point(1, 0, 0);
                nbrPatchNormals[facei] = thisPatchNormals[facei];
            }
            else
            {
                scalar areaDiff = mag(magSf - nbrMagSf)/avSf;

                if (areaDiff > maxAreaDiff)
                {
                    maxAreaDiff = areaDiff;
                    maxAreaFacei = facei;
                }

                if (areaDiff > matchTolerance())
                {
                    FatalErrorInFunction
                        << "face " << facei
                        << " area does not match neighbour by "
                        << 100*areaDiff
                        << "% -- possible face ordering problem." << endl
                        << "patch:" << name()
                        << " my area:" << magSf
                        << " neighbour area:" << nbrMagSf
                        << " matching tolerance:" << matchTolerance()
                         << endl
                        << "Mesh face:" << start()+facei
                        << " fc:" << thisPatchCtrs[facei]
                        << endl
                        << "Neighbour fc:" << nbrPatchCtrs[facei]
                        << endl
                        << "If you are certain your matching is correct"
                        << " you can increase the 'matchTolerance' setting"
                        << " in the patch dictionary in the boundary file."
                        << endl
                        << "Rerun with cyclic debug flag set"
                        << " for more information." << exit(FatalError);
                }
                else
                {
                    thisPatchNormals[facei] = thisPatchAreas[facei] / magSf;
                    nbrPatchNormals[facei] = nbrPatchAreas[facei] / nbrMagSf;
                }
            }
        }


        // Print area match
        if (debug)
        {
            Pout<< "cyclicPolyPatch::calcTransforms :"
                << " patch:" << name()
                << " Max area error:" << 100*maxAreaDiff << "% at face:"
                << maxAreaFacei << " at:" << thisPatchCtrs[maxAreaFacei]
                << " coupled face at:" << nbrPatchCtrs[maxAreaFacei]
                << endl;
        }


        // Calculate transformation

        if (transformType() == ROTATIONAL)
        {
            // Calculate using the given rotation axis and centre. Do not
            // use calculated normals.
            vector n0 = findFaceMaxRadius(thisPatchCtrs);
            vector n1 = -findFaceMaxRadius(nbrPatchCtrs);
            n0 /= mag(n0) + vSmall;
            n1 /= mag(n1) + vSmall;

            if (debug)
            {
                scalar theta = radToDeg(acos(n0 & n1));

                Pout<< "cyclicPolyPatch::calcTransforms :"
                    << " patch:" << name()
                    << " Specified rotation :"
                    << " n0:" << n0 << " n1:" << n1
                    << " swept angle: " << theta << " [deg]"
                    << endl;
            }

            // Extended tensor from two local coordinate systems calculated
            // using normal and rotation axis
            const tensor E0
            (
                rotationAxis_,
                (n0 ^ rotationAxis_),
                n0
            );
            const tensor E1
            (
                rotationAxis_,
                (-n1 ^ rotationAxis_),
                -n1
            );
            const tensor revT(E1.T() & E0);

            transform_ = transformer(revT.T());
        }
        else if (transformType() == TRANSLATIONAL)
        {
            if (debug)
            {
                Pout<< "cyclicPolyPatch::calcTransforms :"
                    << " patch:" << name()
                    << " Specified separation vector : "
                    << separation_ << endl;
            }

            const scalarField thisPatchTols
            (
                matchTolerance()
               *calcFaceTol
                (
                    thisPatch,
                    thisPatch.points(),
                    static_cast<const pointField&>(thisPatchCtrs)
                )
            );

            // Check that separation vectors are same.
            const scalar avgTol = average(thisPatchTols);
            if
            (
                mag(separation_ + nbrPatch().separation_) > avgTol
            )
            {
                WarningInFunction
                    << "Specified separation vector " << separation_
                    << " differs by that of neighbouring patch "
                    << nbrPatch().separation_
                    << " by more than tolerance " << avgTol << endl
                    << "patch:" << name()
                    << " neighbour:" << nbrPatchName()
                    << endl;
            }

            // Set transformation
            transform_ = transformer(separation_);
        }
        else
        {
            const scalarField thisPatchTols
            (
                matchTolerance()
               *calcFaceTol
                (
                    thisPatch,
                    thisPatch.points(),
                    static_cast<const pointField&>(thisPatchCtrs)
                )
            );

            calcTransformTensors
            (
                thisPatchCtrs,
                nbrPatchCtrs,
                thisPatchNormals,
                nbrPatchNormals,
                thisPatchTols,
                matchTolerance(),
                ordering(),
                transformType()
            );
        }
    }
}


void Foam::cyclicPolyPatch::getCentresAndAnchors
(
    const primitivePatch& pp0,
    const primitivePatch& pp1,

    pointField& thisPatchCtrs,
    pointField& nbrPatchCtrs,
    pointField& anchors0,
    scalarField& tols
) const
{
    // Get geometric data on both halves.
    thisPatchCtrs = pp0.faceCentres();
    anchors0 = getAnchorPoints(pp0, pp0.points(), ordering());
    nbrPatchCtrs = pp1.faceCentres();

    if (debug)
    {
        Pout<< "cyclicPolyPatch::getCentresAndAnchors :"
            << " patch:" << name() << nl
            << "thisPatch untransformed faceCentres (avg) : "
            << gAverage(thisPatchCtrs) << nl
            << "nbrPatch untransformed faceCentres (avg) : "
            << gAverage(nbrPatchCtrs) << endl;
    }

    if (thisPatchCtrs.size())
    {
        switch (transformType())
        {
            case ROTATIONAL:
            {
                vector n0 = findFaceMaxRadius(thisPatchCtrs);
                vector n1 = -findFaceMaxRadius(nbrPatchCtrs);
                n0 /= mag(n0) + vSmall;
                n1 /= mag(n1) + vSmall;

                if (debug)
                {
                    scalar theta = radToDeg(acos(n0 & n1));

                    Pout<< "cyclicPolyPatch::getCentresAndAnchors :"
                        << " patch:" << name()
                        << " Specified rotation :"
                        << " n0:" << n0 << " n1:" << n1
                        << " swept angle: " << theta << " [deg]"
                        << endl;
                }

                // Extended tensor from two local coordinate systems calculated
                // using normal and rotation axis
                const tensor E0
                (
                    rotationAxis_,
                    (n0 ^ rotationAxis_),
                    n0
                );
                const tensor E1
                (
                    rotationAxis_,
                    (-n1 ^ rotationAxis_),
                    -n1
                );
                const tensor revT(E1.T() & E0);

                // Rotation
                forAll(thisPatchCtrs, facei)
                {
                    thisPatchCtrs[facei] =
                        Foam::transform
                        (
                            revT,
                            thisPatchCtrs[facei] - rotationCentre_
                        )
                      + rotationCentre_;
                    anchors0[facei] =
                        Foam::transform
                        (
                            revT,
                            anchors0[facei] - rotationCentre_
                        )
                      + rotationCentre_;
                }

                break;
            }
            case TRANSLATIONAL:
            {
                // Transform 0 points.

                if (debug)
                {
                    Pout<< "cyclicPolyPatch::getCentresAndAnchors :"
                        << " patch:" << name()
                        << "Specified translation : " << separation_
                        << endl;
                }

                // Note: getCentresAndAnchors gets called on the slave side
                // so separation is owner-slave points.

                thisPatchCtrs -= separation_;
                anchors0 -= separation_;
                break;
            }
            default:
            {
                // Assumes that cyclic is rotational. This is also the initial
                // condition for patches without faces.

                // Determine the face with max area on both halves. These
                // two faces are used to determine the transformation tensors
                const label max0I = findMaxArea(pp0.points(), pp0);
                const vector n0 = pp0[max0I].normal(pp0.points());

                const label max1I = findMaxArea(pp1.points(), pp1);
                const vector n1 = pp1[max1I].normal(pp1.points());

                if (mag(n0 & n1) < 1-matchTolerance())
                {
                    if (debug)
                    {
                        Pout<< "cyclicPolyPatch::getCentresAndAnchors :"
                            << " patch:" << name()
                            << " Detected rotation :"
                            << " n0:" << n0 << " n1:" << n1 << endl;
                    }

                    // Rotation (around origin)
                    const tensor revT(rotationTensor(n0, -n1));

                    // Rotation
                    forAll(thisPatchCtrs, facei)
                    {
                        thisPatchCtrs[facei] = Foam::transform
                        (
                            revT,
                            thisPatchCtrs[facei]
                        );
                        anchors0[facei] = Foam::transform
                        (
                            revT,
                            anchors0[facei]
                        );
                    }
                }
                else
                {
                    // Parallel translation. Get average of all used points.

                    const point ctr0(sum(pp0.localPoints())/pp0.nPoints());
                    const point ctr1(sum(pp1.localPoints())/pp1.nPoints());

                    if (debug)
                    {
                        Pout<< "cyclicPolyPatch::getCentresAndAnchors :"
                            << " patch:" << name()
                            << " Detected translation :"
                            << " n0:" << n0 << " n1:" << n1
                            << " ctr0:" << ctr0 << " ctr1:" << ctr1 << endl;
                    }

                    thisPatchCtrs += ctr1 - ctr0;
                    anchors0 += ctr1 - ctr0;
                }
                break;
            }
        }
    }

    // Calculate typical distance per face
    tols = matchTolerance()*calcFaceTol(pp1, pp1.points(), nbrPatchCtrs);
}


Foam::vector Foam::cyclicPolyPatch::findFaceMaxRadius
(
    const pointField& faceCentres
) const
{
    // Determine a face furthest away from the axis

    const vectorField n((faceCentres - rotationCentre_) ^ rotationAxis_);

    const scalarField magRadSqr(magSqr(n));

    label facei = findMax(magRadSqr);

    if (debug)
    {
        Info<< "findFaceMaxRadius(const pointField&) : patch: " << name() << nl
            << "    rotFace  = " << facei << nl
            << "    point    = " << faceCentres[facei] << nl
            << "    distance = " << Foam::sqrt(magRadSqr[facei])
            << endl;
    }

    return n[facei];
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::cyclicPolyPatch::cyclicPolyPatch
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
    coupledPolyPatch(name, size, start, index, bm, patchType, ordering),
    nbrPatchName_(word::null),
    nbrPatchID_(-1),
    rotationAxis_(Zero),
    rotationCentre_(Zero),
    coupledPointsPtr_(nullptr),
    coupledEdgesPtr_(nullptr)
{
    // Neighbour patch might not be valid yet so no transformation
    // calculation possible.
}


Foam::cyclicPolyPatch::cyclicPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType,
    const word& nbrPatchName,
    const orderingType ordering
)
:
    coupledPolyPatch(name, size, start, index, bm, patchType, ordering),
    nbrPatchName_(nbrPatchName),
    nbrPatchID_(-1),
    rotationAxis_(Zero),
    rotationCentre_(Zero),
    coupledPointsPtr_(nullptr),
    coupledEdgesPtr_(nullptr)
{
    // Neighbour patch might not be valid yet so no transformation
    // calculation possible.
}


Foam::cyclicPolyPatch::cyclicPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType,
    const orderingType ordering
)
:
    coupledPolyPatch(name, dict, index, bm, patchType, ordering),
    cyclicTransform(dict),
    nbrPatchName_(dict.lookupOrDefault("neighbourPatch", word::null)),
    coupleGroup_(dict),
    nbrPatchID_(-1),
    rotationAxis_(Zero),
    rotationCentre_(Zero),
    coupledPointsPtr_(nullptr),
    coupledEdgesPtr_(nullptr)
{
    if (nbrPatchName_ == word::null && !coupleGroup_.valid())
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "No \"neighbourPatch\" provided." << endl
            << "Is your mesh uptodate with split cyclics?" << endl
            << "Run foamUpgradeCyclics to convert mesh and fields"
            << " to split cyclics." << exit(FatalIOError);
    }

    if (nbrPatchName_ == name)
    {
        FatalIOErrorInFunction(dict)
            << "Neighbour patch name " << nbrPatchName_
            << " cannot be the same as this patch " << name
            << exit(FatalIOError);
    }

    switch (transformType())
    {
        case ROTATIONAL:
        {
            dict.lookup("rotationAxis") >> rotationAxis_;
            dict.lookup("rotationCentre") >> rotationCentre_;

            scalar magRot = mag(rotationAxis_);
            if (magRot < small)
            {
                FatalIOErrorInFunction(dict)
                    << "Illegal rotationAxis " << rotationAxis_ << endl
                    << "Please supply a non-zero vector."
                    << exit(FatalIOError);
            }
            rotationAxis_ /= magRot;

            break;
        }
        case TRANSLATIONAL:
        {
            dict.lookup("separation") >> separation_;
            break;
        }
        default:
        {
            // no additional info required
        }
    }

    // Neighbour patch might not be valid yet so no transformation
    // calculation possible.
}


Foam::cyclicPolyPatch::cyclicPolyPatch
(
    const cyclicPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(pp, bm),
    cyclicTransform(pp),
    nbrPatchName_(pp.nbrPatchName_),
    coupleGroup_(pp.coupleGroup_),
    nbrPatchID_(-1),
    rotationAxis_(pp.rotationAxis_),
    rotationCentre_(pp.rotationCentre_),
    coupledPointsPtr_(nullptr),
    coupledEdgesPtr_(nullptr)
{
    // Neighbour patch might not be valid yet so no transformation
    // calculation possible.
}


Foam::cyclicPolyPatch::cyclicPolyPatch
(
    const cyclicPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart,
    const word& neiName
)
:
    coupledPolyPatch(pp, bm, index, newSize, newStart),
    cyclicTransform(pp),
    nbrPatchName_(neiName),
    coupleGroup_(pp.coupleGroup_),
    nbrPatchID_(-1),
    rotationAxis_(pp.rotationAxis_),
    rotationCentre_(pp.rotationCentre_),
    coupledPointsPtr_(nullptr),
    coupledEdgesPtr_(nullptr)
{
    if (neiName == name())
    {
        FatalErrorInFunction
            << "Neighbour patch name " << neiName
            << " cannot be the same as this patch " << name()
            << exit(FatalError);
    }

    // Neighbour patch might not be valid yet so no transformation
    // calculation possible.
}


Foam::cyclicPolyPatch::cyclicPolyPatch
(
    const cyclicPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    coupledPolyPatch(pp, bm, index, mapAddressing, newStart),
    cyclicTransform(pp),
    nbrPatchName_(pp.nbrPatchName_),
    coupleGroup_(pp.coupleGroup_),
    nbrPatchID_(-1),
    rotationAxis_(pp.rotationAxis_),
    rotationCentre_(pp.rotationCentre_),
    coupledPointsPtr_(nullptr),
    coupledEdgesPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cyclicPolyPatch::~cyclicPolyPatch()
{
    deleteDemandDrivenData(coupledPointsPtr_);
    deleteDemandDrivenData(coupledEdgesPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::word& Foam::cyclicPolyPatch::nbrPatchName() const
{
    if (nbrPatchName_.empty())
    {
        // Try and use patchGroup to find samplePatch and sampleRegion
        label patchID = coupleGroup_.findOtherPatchID(*this);

        nbrPatchName_ = boundaryMesh()[patchID].name();
    }
    return nbrPatchName_;
}


Foam::label Foam::cyclicPolyPatch::nbrPatchID() const
{
    if (nbrPatchID_ == -1)
    {
        nbrPatchID_ = this->boundaryMesh().findPatchID(nbrPatchName());

        if (nbrPatchID_ == -1)
        {
            FatalErrorInFunction
                << "Illegal neighbourPatch name " << nbrPatchName()
                << endl << "Valid patch names are "
                << this->boundaryMesh().names()
                << exit(FatalError);
        }

        // Check that it is a cyclic
        const cyclicPolyPatch& nbrPatch = refCast<const cyclicPolyPatch>
        (
            this->boundaryMesh()[nbrPatchID_]
        );

        if (nbrPatch.nbrPatchName() != name())
        {
            WarningInFunction
                << "Patch " << name()
                << " specifies neighbour patch " << nbrPatchName()
                << endl << " but that in return specifies "
                << nbrPatch.nbrPatchName()
                << endl;
        }
    }
    return nbrPatchID_;
}


void Foam::cyclicPolyPatch::transformPosition(pointField& l) const
{
    if (transform().rotates())
    {
        if (transformType() == ROTATIONAL)
        {
            l =
                Foam::transform(transform().R(), l-rotationCentre_)
              + rotationCentre_;
        }
        else
        {
            l = Foam::transform(transform().R(), l);
        }
    }
    else if (transform().translates())
    {
        // transformPosition gets called on the receiving side,
        // separation gets calculated on the sending side so subtract.
        l -= transform().t();
    }
}


void Foam::cyclicPolyPatch::transformPosition(point& l, const label facei) const
{
    if (transform().rotates())
    {
        if (transformType() == ROTATIONAL)
        {
            l = Foam::transform(transform().R(), l - rotationCentre_)
              + rotationCentre_;
        }
        else
        {
            l = Foam::transform(transform().R(), l);
        }
    }
    else if (transform().translates())
    {
        l -= transform().t();
    }
}


void Foam::cyclicPolyPatch::initGeometry(PstreamBuffers& pBufs)
{
    polyPatch::initGeometry(pBufs);
}


void Foam::cyclicPolyPatch::initGeometry
(
    const primitivePatch& referPatch,
    pointField& nbrCtrs,
    vectorField& nbrAreas,
    pointField& nbrCc
)
{}


void Foam::cyclicPolyPatch::calcGeometry
(
    const primitivePatch& referPatch,
    const pointField& thisCtrs,
    const vectorField& thisAreas,
    const pointField& thisCc,
    const pointField& nbrCtrs,
    const vectorField& nbrAreas,
    const pointField& nbrCc
)
{
    calcTransforms
    (
        referPatch,
        thisCtrs,
        thisAreas,
        nbrCtrs,
        nbrAreas
    );
}


void Foam::cyclicPolyPatch::calcGeometry(PstreamBuffers& pBufs)
{
    calcGeometry
    (
        *this,
        faceCentres(),
        faceAreas(),
        faceCellCentres(),
        nbrPatch().faceCentres(),
        nbrPatch().faceAreas(),
        nbrPatch().faceCellCentres()
    );
}


void Foam::cyclicPolyPatch::initMovePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    polyPatch::initMovePoints(pBufs, p);
}


void Foam::cyclicPolyPatch::movePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    polyPatch::movePoints(pBufs, p);
    calcTransforms();
}


void Foam::cyclicPolyPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    polyPatch::initUpdateMesh(pBufs);
}


void Foam::cyclicPolyPatch::updateMesh(PstreamBuffers& pBufs)
{
    polyPatch::updateMesh(pBufs);
    deleteDemandDrivenData(coupledPointsPtr_);
    deleteDemandDrivenData(coupledEdgesPtr_);
}


const Foam::edgeList& Foam::cyclicPolyPatch::coupledPoints() const
{
    if (!coupledPointsPtr_)
    {
        const faceList& nbrLocalFaces = nbrPatch().localFaces();
        const labelList& nbrMeshPoints = nbrPatch().meshPoints();

        // Now all we know is that relative face index in *this is same
        // as coupled face in nbrPatch and also that the 0th vertex
        // corresponds.

        // From local point to nbrPatch or -1.
        labelList coupledPoint(nPoints(), -1);

        forAll(*this, patchFacei)
        {
            const face& fA = localFaces()[patchFacei];
            const face& fB = nbrLocalFaces[patchFacei];

            forAll(fA, indexA)
            {
                label patchPointA = fA[indexA];

                if (coupledPoint[patchPointA] == -1)
                {
                    label indexB = (fB.size() - indexA) % fB.size();

                    // Filter out points on wedge axis
                    if (meshPoints()[patchPointA] != nbrMeshPoints[fB[indexB]])
                    {
                        coupledPoint[patchPointA] = fB[indexB];
                    }
                }
            }
        }

        coupledPointsPtr_ = new edgeList(nPoints());
        edgeList& connected = *coupledPointsPtr_;

        // Extract coupled points.
        label connectedI = 0;

        forAll(coupledPoint, i)
        {
            if (coupledPoint[i] != -1)
            {
                connected[connectedI++] = edge(i, coupledPoint[i]);
            }
        }

        connected.setSize(connectedI);

        if (debug)
        {
            OFstream str
            (
                boundaryMesh().mesh().time().path()
               /name() + "_coupledPoints.obj"
            );
            label vertI = 0;

            Pout<< "Writing file " << str.name() << " with coordinates of "
                << "coupled points" << endl;

            forAll(connected, i)
            {
                const point& a = points()[meshPoints()[connected[i][0]]];
                const point& b = points()[nbrMeshPoints[connected[i][1]]];

                str<< "v " << a.x() << ' ' << a.y() << ' ' << a.z() << nl;
                str<< "v " << b.x() << ' ' << b.y() << ' ' << b.z() << nl;
                vertI += 2;

                str<< "l " << vertI-1 << ' ' << vertI << nl;
            }
        }
    }
    return *coupledPointsPtr_;
}


const Foam::edgeList& Foam::cyclicPolyPatch::coupledEdges() const
{
    if (!coupledEdgesPtr_)
    {
        const edgeList& pointCouples = coupledPoints();

        // Build map from points on *this (A) to points on neighbourpatch (B)
        Map<label> aToB(2*pointCouples.size());

        forAll(pointCouples, i)
        {
            const edge& e = pointCouples[i];

            aToB.insert(e[0], e[1]);
        }

        // Map from edge on A to points (in B indices)
        EdgeMap<label> edgeMap(nEdges());

        forAll(*this, patchFacei)
        {
            const labelList& fEdges = faceEdges()[patchFacei];

            forAll(fEdges, i)
            {
                label edgeI = fEdges[i];

                const edge& e = edges()[edgeI];

                // Convert edge end points to corresponding points on B side.
                Map<label>::const_iterator fnd0 = aToB.find(e[0]);
                if (fnd0 != aToB.end())
                {
                    Map<label>::const_iterator fnd1 = aToB.find(e[1]);
                    if (fnd1 != aToB.end())
                    {
                        edgeMap.insert(edge(fnd0(), fnd1()), edgeI);
                    }
                }
            }
        }

        // Use the edgeMap to get the edges on the B side.

        const cyclicPolyPatch& nbrPatch = this->nbrPatch();
        const labelList& nbrMp = nbrPatch.meshPoints();
        const labelList& mp = meshPoints();



        coupledEdgesPtr_ = new edgeList(edgeMap.size());
        edgeList& coupledEdges = *coupledEdgesPtr_;
        label coupleI = 0;

        forAll(nbrPatch, patchFacei)
        {
            const labelList& fEdges = nbrPatch.faceEdges()[patchFacei];

            forAll(fEdges, i)
            {
                label edgeI = fEdges[i];

                const edge& e = nbrPatch.edges()[edgeI];

                // Look up A edge from HashTable.
                EdgeMap<label>::iterator iter = edgeMap.find(e);

                if (iter != edgeMap.end())
                {
                    label edgeA = iter();
                    const edge& eA = edges()[edgeA];

                    // Store correspondence. Filter out edges on wedge axis.
                    if
                    (
                        edge(mp[eA[0]], mp[eA[1]])
                     != edge(nbrMp[e[0]], nbrMp[e[1]])
                    )
                    {
                        coupledEdges[coupleI++] = edge(edgeA, edgeI);
                    }

                    // Remove so we build unique list
                    edgeMap.erase(iter);
                }
            }
        }
        coupledEdges.setSize(coupleI);


        // Some checks

        forAll(coupledEdges, i)
        {
            const edge& e = coupledEdges[i];

            if (e[0] < 0 || e[1] < 0)
            {
                FatalErrorInFunction
                    << "Problem : at position " << i
                    << " illegal couple:" << e
                    << abort(FatalError);
            }
        }

        if (debug)
        {
            OFstream str
            (
                boundaryMesh().mesh().time().path()
               /name() + "_coupledEdges.obj"
            );
            label vertI = 0;

            Pout<< "Writing file " << str.name() << " with centres of "
                << "coupled edges" << endl;

            forAll(coupledEdges, i)
            {
                const edge& e = coupledEdges[i];

                const point& a = edges()[e[0]].centre(localPoints());
                const point& b = nbrPatch.edges()[e[1]].centre
                (
                    nbrPatch.localPoints()
                );

                str<< "v " << a.x() << ' ' << a.y() << ' ' << a.z() << nl;
                str<< "v " << b.x() << ' ' << b.y() << ' ' << b.z() << nl;
                vertI += 2;

                str<< "l " << vertI-1 << ' ' << vertI << nl;
            }
        }
    }
    return *coupledEdgesPtr_;
}


void Foam::cyclicPolyPatch::initOrder
(
    PstreamBuffers&,
    const primitivePatch& pp
) const
{
    if (owner())
    {
        // Save patch for use in non-owner side ordering. Equivalent to
        // processorPolyPatch using OPstream.
        ownerPatchPtr_.reset
        (
            new primitivePatch
            (
                pp,
                pp.points()
            )
        );
    }
}


bool Foam::cyclicPolyPatch::order
(
    PstreamBuffers& pBufs,
    const primitivePatch& pp,
    labelList& faceMap,
    labelList& rotation
) const
{
    if (debug)
    {
        Pout<< "order : of " << pp.size()
            << " faces of patch:" << name()
            << " neighbour:" << nbrPatchName()
            << endl;
    }
    faceMap.setSize(pp.size());
    faceMap = -1;

    rotation.setSize(pp.size());
    rotation = 0;

    if (ordering() == NOORDERING)
    {
        // No faces, nothing to change.
        return false;
    }

    if (owner())
    {
        // Do nothing (i.e. identical mapping, zero rotation).
        // See comment at top.
        forAll(faceMap, patchFacei)
        {
            faceMap[patchFacei] = patchFacei;
        }

        return false;
    }
    else
    {
        // Get stored geometry from initOrder invocation of owner.
        const primitivePatch& pp0 = nbrPatch().ownerPatchPtr_();

        // Get geometric quantities
        pointField thisPatchCtrs, nbrPatchCtrs, anchors0;
        scalarField tols;
        getCentresAndAnchors
        (
            pp0,
            pp,

            thisPatchCtrs,
            nbrPatchCtrs,
            anchors0,
            tols
        );

        if (debug)
        {
            Pout<< "thisPatch transformed faceCentres (avg)   : "
                << gAverage(thisPatchCtrs) << nl
                << "nbrPatch untransformed faceCentres (avg) : "
                << gAverage(nbrPatchCtrs) << endl;
        }

        // Geometric match of face centre vectors
        bool matchedAll = matchPoints
        (
            nbrPatchCtrs,
            thisPatchCtrs,
            tols,
            true,
            faceMap
        );

        if (!matchedAll || debug)
        {
            // Dump halves
            fileName nm0
            (
                boundaryMesh().mesh().time().path()
               /nbrPatch().name()+"_faces.obj"
            );
            Pout<< "cyclicPolyPatch::order : Writing neighbour"
                << " faces to OBJ file " << nm0 << endl;
            writeOBJ(nm0, pp0, pp0.points());

            fileName nm1
            (
                boundaryMesh().mesh().time().path()
               /name()+"_faces.obj"
            );
            Pout<< "cyclicPolyPatch::order : Writing my"
                << " faces to OBJ file " << nm1 << endl;
            writeOBJ(nm1, pp, pp.points());

            OFstream ccStr
            (
                boundaryMesh().mesh().time().path()
               /name() + "_faceCentres.obj"
            );
            Pout<< "cyclicPolyPatch::order : "
                << "Dumping currently found cyclic match as lines between"
                << " corresponding face centres to file " << ccStr.name()
                << endl;

            // Recalculate untransformed face centres
            // pointField rawthisPatchCtrs =
            //    calcFaceCentres(thisPatchFaces, pp.points());
            label vertI = 0;

            forAll(nbrPatchCtrs, i)
            {
                if (faceMap[i] != -1)
                {
                    // Write edge between c1 and c0
                    const point& c0 = thisPatchCtrs[faceMap[i]];
                    const point& c1 = nbrPatchCtrs[i];
                    writeOBJ(ccStr, c0, c1, vertI);
                }
            }
        }

        if (!matchedAll)
        {
            SeriousErrorInFunction
                << "Patch:" << name() << " : "
                << "Cannot match vectors to faces on both sides of patch"
                << endl
                << "    Perhaps your faces do not match?"
                << " The obj files written contain the current match." << endl
                << "    Continuing with incorrect face ordering from now on!"
                << endl;

                return false;
        }


        // Set rotation.
        forAll(faceMap, oldFacei)
        {
            // The face f will be at newFacei (after morphing) and we want its
            // anchorPoint (= f[0]) to align with the anchorpoint for the
            // corresponding face on the other side.

            label newFacei = faceMap[oldFacei];

            const point& wantedAnchor = anchors0[newFacei];

            rotation[newFacei] = getRotation
            (
                pp.points(),
                pp[oldFacei],
                wantedAnchor,
                tols[oldFacei]
            );

            if (rotation[newFacei] == -1)
            {
                SeriousErrorInFunction
                    << "in patch " << name()
                    << " : "
                    << "Cannot find point on face " << pp[oldFacei]
                    << " with vertices "
                    << IndirectList<point>(pp.points(), pp[oldFacei])()
                    << " that matches point " << wantedAnchor
                    << " when matching the halves of processor patch " << name()
                    << "Continuing with incorrect face ordering from now on!"
                    << endl;

                return false;
            }
        }

        ownerPatchPtr_.clear();

        // Return false if no change necessary, true otherwise.

        forAll(faceMap, facei)
        {
            if (faceMap[facei] != facei || rotation[facei] != 0)
            {
                return true;
            }
        }

        return false;
    }
}


void Foam::cyclicPolyPatch::write(Ostream& os) const
{
    coupledPolyPatch::write(os);

    if (!nbrPatchName_.empty())
    {
        writeEntry(os, "neighbourPatch", nbrPatchName_);
    }

    coupleGroup_.write(os);

    cyclicTransform::write(os);

    switch (transformType())
    {
        case ROTATIONAL:
        {
            writeEntry(os, "rotationAxis", rotationAxis_);
            writeEntry(os, "rotationCentre", rotationCentre_);
            break;
        }
        case TRANSLATIONAL:
        {
            writeEntry(os, "separation", separation_);
            break;
        }
        default:
        {
            // no additional info to write
        }
    }
}


// ************************************************************************* //
