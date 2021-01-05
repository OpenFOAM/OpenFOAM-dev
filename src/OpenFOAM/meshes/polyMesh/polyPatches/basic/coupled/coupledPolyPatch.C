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

#include "coupledPolyPatch.H"
#include "ListOps.H"
#include "transform.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(coupledPolyPatch, 0);

    const scalar coupledPolyPatch::defaultMatchTol_ = 1e-4;
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

bool Foam::coupledPolyPatch::walk
(
    const primitivePatch& pp,
    const bool direction,
    const label seedFacei,
    const label seedFacePointi,
    labelList& faceMap,
    labelList& facePointMap,
    label& mapFacei,
    autoPtr<labelListList>& walks
) const
{
    // Initialisation
    label facei = seedFacei, facePointi = seedFacePointi;
    faceMap[facei] = mapFacei;
    facePointMap[facei] = facePointi;
    bool changed = facei != mapFacei || facePointi != 0;
    ++ mapFacei;
    if (walks.valid())
    {
        walks->append(labelList(1, facei));
    }

    // Walk the patch until we get back to the seed face and point
    do
    {
        // Get the next point around the face
        const label facePointi1 =
            direction
          ? pp[facei].fcIndex(facePointi)
          : pp[facei].rcIndex(facePointi);

        // Get the current edge within the face
        const label faceEdgei =
            direction ? facePointi : pp[facei].rcIndex(facePointi);
        const label edgei = pp.faceEdges()[facei][faceEdgei];

        // If the number of faces connected to this edge is not 2, then this
        // edge is non-manifold and is considered a boundary to the walk. So,
        // the walk moves on to the next point around the current face.
        if (pp.edgeFaces()[edgei].size() != 2)
        {
            facePointi = facePointi1;
            continue;
        }

        // Get the connected face and the corresponding point index within
        const label facej =
            pp.edgeFaces()[edgei][pp.edgeFaces()[edgei][0] == facei];
        const label facePointj = pp[facej].which(pp[facei][facePointi]);

        // Get the corresponding next point within the connected face
        const label facePointj1 =
            direction
          ? pp[facej].rcIndex(facePointj)
          : pp[facej].fcIndex(facePointj);

        // If the next points are not the same then that indicates that the
        // faces are numbered in opposite directions. This means that the faces
        // are not actually connected. There should really be two edges, each
        // connected to just one of the faces. This edge should therefore be
        // considered a boundary to the walk.
        if (pp[facei][facePointi1] != pp[facej][facePointj1])
        {
            facePointi = facePointi1;
            continue;
        }

        // It has been determined that the current edge *can* be crossed. Now
        // test whether of not the walk *should* cross this edge...
        if (faceMap[facej] == -1)
        {
            // The connected face has not been visited, so walk into it and
            // set its ordering in the map, its point and visited status
            facei = facej;
            facePointi = facePointj;
            faceMap[facei] = mapFacei;
            facePointMap[facei] = facePointi;
            changed = changed || facei != mapFacei || facePointi != 0;
            ++ mapFacei;
        }
        else if (facePointMap[facei] != facePointi1 || facei == seedFacei)
        {
            // The connected face has been visited, but there are more
            // edges to consider on the current face, so move to the next
            // face point
            facePointi = facePointi1;
        }
        else
        {
            // The connected face has been visited, and there are no more
            // edges to consider on the current face, so backtrack to the
            // previous face in the walk
            facei = facej;
            facePointi = facePointj;
        }

        // Add to the walk, if that information is being stored
        if (walks.valid() && walks->last().last() != facei)
        {
            walks->last().append(facei);
        }
    }
    while (facei != seedFacei || facePointi != seedFacePointi);

    return changed;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::coupledPolyPatch::writeOBJ
(
    const fileName& name,
    const primitivePatch& pp
)
{
    OFstream os(name);

    forAll(pp.localPoints(), pointi)
    {
        const point& p = pp.localPoints()[pointi];
        os << "v " << p.x() << ' ' << p.y() << ' ' << p.z() << endl;
    }

    forAll(pp.localFaces(), facei)
    {
        const face& f = pp.localFaces()[facei];
        os << 'f';
        forAll(f, fi)
        {
            os << ' ' << f[fi] + 1;
        }
        os << endl;
    }
}

void Foam::coupledPolyPatch::writeOBJ
(
    const fileName& name,
    const pointField& points0,
    const pointField& points1
)
{
    OFstream os(name);

    forAll(points0, pointi)
    {
        const point& p0 = points0[pointi];
        const point& p1 = points1[pointi];
        os  << "v " << p0.x() << ' ' << p0.y() << ' ' << p0.z() << endl
            << "v " << p1.x() << ' ' << p1.y() << ' ' << p1.z() << endl
            << "l " << 2*pointi << ' ' << 2*pointi + 1 << endl;
    }
}


void Foam::coupledPolyPatch::writeOBJ
(
    const fileName& name,
    const pointField& points,
    const labelListList& paths
)
{
    OFstream os(name);

    forAll(points, pointi)
    {
        const point& c = points[pointi];
        os << "v " << c.x() << ' '<< c.y() << ' ' << c.z() << endl;
    }

    forAll(paths, pathi)
    {
        for (label pathj = 0; pathj < paths[pathi].size() - 1; ++ pathj)
        {
            os  << "l " << paths[pathi][pathj] + 1 << ' '
                << paths[pathi][pathj + 1] + 1 << endl;
        }
    }
}


void Foam::coupledPolyPatch::initOrder
(
    ownToNbrOrderData& ownToNbr,
    autoPtr<ownToNbrDebugOrderData>& ownToNbrDebugPtr,
    const primitivePatch& pp
) const
{
    if (owner())
    {
        // Generate the connected regions
        label nRegions = 0;
        labelList faceRegionis(pp.size(), -1);

        label seedFacei = 0;

        labelList faceMap(pp.size(), -1);
        labelList facePointMap(pp.size(), -1);
        label mapFacei = 0;
        autoPtr<labelListList> walks(nullptr);

        while (mapFacei < pp.size())
        {
            walk
            (
                pp,
                owner(),
                seedFacei,
                0,
                faceMap,
                facePointMap,
                mapFacei,
                walks
            );

            forAll(pp, facei)
            {
                if (faceRegionis[facei] == -1 && faceMap[facei] != -1)
                {
                    faceRegionis[facei] = nRegions;
                }
            }

            ++ nRegions;

            forAll(pp, facei)
            {
                if (faceMap[facei] == -1)
                {
                    seedFacei = facei;
                    break;
                }
            }
        }

        // Generate the face tolerances
        //
        // !!! It is possble that a different metric would be more appropriate
        // for this method than the tolerance that was developed when all faces
        // were being geometrically compared
        //
        const scalarField tols(calcFaceTol(pp, pp.points(), pp.faceCentres()));

        // Get the face with the largest tolerance in each region as the seed
        // and store its index in the (self) ordering data
        ownToOwnOrderDataPtr_ = new ownToOwnOrderData();
        ownToOwnOrderDataPtr_->seedFaceis = labelList(nRegions, -1);
        scalarList maxTols(nRegions, -vGreat);
        forAll(pp, facei)
        {
            const label regioni = faceRegionis[facei];

            if (tols[facei] > maxTols[regioni])
            {
                ownToOwnOrderDataPtr_->seedFaceis[regioni] = facei;
                maxTols[regioni] = tols[facei];
            }
        }

        // Get the points of each seed face and store them in the neighbour
        // ordering data
        ownToNbr.seedFacePoints = List<pointField>(nRegions, pointField());
        forAll(ownToOwnOrderDataPtr_->seedFaceis, regioni)
        {
            const face& f = pp[ownToOwnOrderDataPtr_->seedFaceis[regioni]];
            ownToNbr.seedFacePoints[regioni] = f.points(pp.points());
        }

        // Get debug data
        if (ownToNbrDebugPtr.valid())
        {
            ownToNbrDebugPtr = new ownToNbrDebugOrderData();
            ownToNbrDebugPtr->nFaces = pp.size();
            ownToNbrDebugPtr->nPoints = pp.nPoints();
            ownToNbrDebugPtr->nEdges = pp.nEdges();
            ownToNbrDebugPtr->nInternalEdges = pp.nInternalEdges();
        }
    }
}


bool Foam::coupledPolyPatch::order
(
    const ownToNbrOrderData& ownToNbr,
    const autoPtr<ownToNbrDebugOrderData>& ownToNbrDebugPtr,
    const primitivePatch& pp,
    labelList& faceMap,
    labelList& rotation
) const
{
    // Determine the seed faces and face points
    labelList seedFaceis, seedFacePointis;
    if (owner())
    {
        seedFaceis = ownToOwnOrderDataPtr_->seedFaceis;
        ownToOwnOrderDataPtr_.clear();

        seedFacePointis.resize(seedFaceis.size(), 0);
    }
    else
    {
        const List<pointField> ownerSeedFacePoints(ownToNbr.seedFacePoints);

        seedFaceis.resize(ownerSeedFacePoints.size());
        seedFacePointis.resize(ownerSeedFacePoints.size());

        // Check the element counts
        if (ownToNbrDebugPtr.valid())
        {
            const label ownerNFaces = ownToNbrDebugPtr->nFaces;
            const label ownerNPoints = ownToNbrDebugPtr->nPoints;
            const label ownerNEdges = ownToNbrDebugPtr->nEdges;
            const label ownerNInternalEdges = ownToNbrDebugPtr->nInternalEdges;
            if (pp.size() != ownerNFaces)
            {
                SeriousErrorInFunction<< "The patch " << name() << " has "
                    << pp.size() << " faces whilst it's neighbour has "
                    << ownerNFaces << endl;
            }
            if (pp.nPoints() != ownerNPoints)
            {
                SeriousErrorInFunction<< "The patch " << name() << " has "
                    << pp.nPoints() << " points whilst it's neighbour has "
                    << ownerNPoints << endl;
            }
            if (pp.nEdges() != ownerNEdges)
            {
                SeriousErrorInFunction<< "The patch " << name() << " has "
                    << pp.nEdges() << " edges whilst it's neighbour has "
                    << ownerNEdges << endl;
            }
            if (pp.nInternalEdges() != ownerNInternalEdges)
            {
                SeriousErrorInFunction<< "The patch " << name() << " has "
                    << pp.nInternalEdges() << " internal edges whilst it's "
                    << "neighbour has " << ownerNInternalEdges << endl;
            }
        }

        // Do geometric testing to determine the faces that match those sent
        // over from the opposite patch
        forAll(ownerSeedFacePoints, regioni)
        {
            const pointField& ownerFacePts = ownerSeedFacePoints[regioni];

            // The seed face and face-point are the ones which give the smallest
            // total displacement between all corresponding points. Note that
            // owner and neighbour point order is reversed.
            scalar minSumSqrDisplacement = vGreat;
            forAll(pp, facei)
            {
                const pointField facePts = pp[facei].points(pp.points());

                if (facePts.size() != ownerFacePts.size()) continue;

                forAll(facePts, facePointi)
                {
                    const scalar sumSqrDisplacement =
                        sum
                        (
                            magSqr
                            (
                                rotateList(reverseList(facePts), facePointi + 1)
                              - ownerFacePts
                            )
                        );
                    if (sumSqrDisplacement < minSumSqrDisplacement)
                    {
                        seedFaceis[regioni] = facei;
                        seedFacePointis[regioni] = facePointi;
                        minSumSqrDisplacement = sumSqrDisplacement;
                    }
                }
            }

            // Check and report if the min displacement is large
            const scalar seedFaceTol =
                calcFaceTol
                (
                    faceList(1, pp[seedFaceis[regioni]]),
                    pp.points(),
                    pointField(1, pp.faceCentres()[seedFaceis[regioni]])
                ).first();
            if (minSumSqrDisplacement > seedFaceTol)
            {
                FatalErrorInFunction
                    << "The root-sum-square displacement between the points of "
                    << "the seed face and the best matching face (#"
                    << seedFaceis[regioni] << ") on patch " << name()
                    << " is " << sqrt(minSumSqrDisplacement) << "."
                    << nl
                    << "This is greater than the match tolerance of "
                    << seedFaceTol << " for this face."
                    << nl
                    << "Check that the patches are conformal and that any "
                    << "transformations defined between them are correct"
                    << nl
                    << "It might be possible to fix this problem by increasing "
                    << "the \"matchTolerance\" setting for this patch in the "
                    << "boundary file."
                    << nl
                    << "Re-run with the \"coupled\" debug flag set for more "
                    << "information."
                    << exit(FatalError);
            }
        }
    }

    // Walk the patch from the seeds
    bool changed = false;
    faceMap = -1;
    labelList facePointMap(pp.size(), -1);
    label mapFacei = 0;
    autoPtr<labelListList> walks(debug ? new labelListList() : nullptr);
    forAll(seedFaceis, regioni)
    {
        changed =
            walk
            (
                pp,
                owner(),
                seedFaceis[regioni],
                seedFacePointis[regioni],
                faceMap,
                facePointMap,
                mapFacei,
                walks
            )
         || changed;
    }

    // Write out the patch
    if (debug)
    {
        Pout<< "Writing patch " << name() << " to " << name() + ".obj" << endl;
        writeOBJ(name() + ".obj", pp);
    }

    // Write out the walk
    if (debug)
    {
        Pout<< "Writing patch " << name() << " walks to "
            << name() + "Walk.obj" << endl;
        writeOBJ(name() + "Walk.obj", pp.faceCentres(), walks());
    }

    // Check that all faces have been visited exactly once
    bool badWalk = mapFacei != pp.size();
    forAll(pp, facei)
    {
        badWalk = badWalk || facePointMap[facei] == -1;
    }
    if (badWalk)
    {
        FatalErrorInFunction
            << "The ordering walk did not hit every face exactly once"
            << exit(FatalError);
    }

    // Construct the rotations from the face point map
    forAll(pp, facei)
    {
        rotation[facei] =
            (pp[facei].size() - facePointMap[facei]) % pp[facei].size();
    }

    // Map the rotations
    //
    // !!! The rotation list appears to be indexed by the new face label,
    // rather than the old one. For sanity's sake the ordering code above
    // indexes everything consistently with the old face label. This means the
    // rotations need mapping to the new indices.
    //
    UIndirectList<label>(rotation, faceMap) = labelList(rotation);

    return changed;
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::coupledPolyPatch::coupledPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    polyPatch(name, size, start, index, bm, patchType),
    matchTolerance_(defaultMatchTol_),
    ownToOwnOrderDataPtr_(nullptr)
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
    ownToOwnOrderDataPtr_(nullptr)
{}


Foam::coupledPolyPatch::coupledPolyPatch
(
    const coupledPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    polyPatch(pp, bm),
    matchTolerance_(pp.matchTolerance_),
    ownToOwnOrderDataPtr_(nullptr)
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
    ownToOwnOrderDataPtr_(nullptr)
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
    ownToOwnOrderDataPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coupledPolyPatch::~coupledPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

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


void Foam::coupledPolyPatch::write(Ostream& os) const
{
    polyPatch::write(os);
    writeEntry(os, "matchTolerance", matchTolerance_);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>
(
    Istream& is,
    coupledPolyPatch::ownToNbrOrderData& ownToNbr
)
{
    is >> ownToNbr.seedFacePoints;
    return is;
}


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const coupledPolyPatch::ownToNbrOrderData& ownToNbr
)
{
    os << ownToNbr.seedFacePoints;
    return os;
}


Foam::Istream& Foam::operator>>
(
    Istream& is,
    coupledPolyPatch::ownToNbrDebugOrderData& ownToNbrDebug
)
{
    is  >> ownToNbrDebug.nFaces
        >> ownToNbrDebug.nPoints
        >> ownToNbrDebug.nEdges
        >> ownToNbrDebug.nInternalEdges;
    return is;
}


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const coupledPolyPatch::ownToNbrDebugOrderData& ownToNbrDebug
)
{
    os  << ownToNbrDebug.nFaces
        << ownToNbrDebug.nPoints
        << ownToNbrDebug.nEdges
        << ownToNbrDebug.nInternalEdges;
    return os;
}


// ************************************************************************* //
