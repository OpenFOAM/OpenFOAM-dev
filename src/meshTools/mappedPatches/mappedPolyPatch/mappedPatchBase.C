/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "mappedPatchBase.H"
#include "meshSearchMeshObject.H"
#include "meshTools.H"
#include "OBJstream.H"
#include "treeDataFace.H"
#include "treeDataCell.H"
#include "indexedOctree.H"
#include "polyMesh.H"
#include "polyPatch.H"
#include "SubField.H"
#include "Time.H"
#include "distributionMap.H"
#include "triPointRef.H"
#include "RemoteData.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mappedPatchBase, 0);

    template<>
    const char* Foam::NamedEnum
    <
        Foam::mappedPatchBase::sampleMode,
        4
    >::names[] =
    {
        "nearestCell",
        "nearestPatchFace",
        "nearestPatchFaceAMI",
        "nearestFace"
    };

    template<>
    const char* Foam::NamedEnum
    <
        Foam::mappedPatchBase::offsetMode,
        3
    >::names[] =
    {
        "none",
        "normal",
        "direction"
    };
}


const Foam::NamedEnum<Foam::mappedPatchBase::sampleMode, 4>
    Foam::mappedPatchBase::sampleModeNames_;

const Foam::NamedEnum<Foam::mappedPatchBase::offsetMode, 3>
    Foam::mappedPatchBase::offsetModeNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::mappedPatchBase::offsetMode Foam::mappedPatchBase::readOffsetMode
(
    const dictionary& dict
) const
{
    if (dict.found("offsetMode"))
    {
        return offsetModeNames_.read(dict.lookup("offsetMode"));
    }
    else
    {
        const bool haveDistance = dict.found("distance");
        const bool haveOffset = dict.found("offset");

        if (haveDistance && haveOffset)
        {
            // Error. Demand "offsetMode" setting to disambiguate.
            return offsetModeNames_.read(dict.lookup("offsetMode"));
        }
        else if (haveDistance)
        {
            return NORMAL;
        }
        else if (haveOffset)
        {
            return DIRECTION;
        }
        else
        {
            return NONE;
        }
    }
}


void Foam::mappedPatchBase::findSamples
(
    const globalIndex& patchGlobalIndex,
    labelList& sampleGlobalPatchFaces,
    labelList& sampleIndices
) const
{
    // Lookup the correct region
    const polyMesh& mesh = sampleMesh();

    // Gather the sample points into a single globally indexed list
    List<point> allPoints(patchGlobalIndex.size());
    {
        List<pointField> procSamplePoints(Pstream::nProcs());
        procSamplePoints[Pstream::myProcNo()] = this->samplePoints();
        Pstream::gatherList(procSamplePoints);
        Pstream::scatterList(procSamplePoints);

        forAll(procSamplePoints, proci)
        {
            forAll(procSamplePoints[proci], procSamplei)
            {
                allPoints
                [
                    patchGlobalIndex.toGlobal(proci, procSamplei)
                ] = procSamplePoints[proci][procSamplei];
            }
        }
    }

    // Nearest info
    List<RemoteData<scalar>> allNearest(patchGlobalIndex.size());

    switch (mode_)
    {
        case NEARESTCELL:
        {
            if (samplePatch_.size() && samplePatch_ != "none")
            {
                FatalErrorInFunction
                    << "No need to supply a patch name when in "
                    << sampleModeNames_[mode_] << " mode." << exit(FatalError);
            }

            //- Note: face-diagonal decomposition
            const indexedOctree<Foam::treeDataCell>& tree = mesh.cellTree();

            forAll(allPoints, alli)
            {
                const point& p = allPoints[alli];

                const label celli = tree.findInside(p);

                if (celli != -1)
                {
                    const point& cc = mesh.cellCentres()[celli];

                    allNearest[alli].proci = Pstream::myProcNo();
                    allNearest[alli].elementi = celli;
                    allNearest[alli].data = magSqr(cc - p);
                }
            }
            break;
        }

        case NEARESTPATCHFACE:
        {
            const polyPatch& pp = samplePolyPatch();

            if (pp.empty())
            {
                forAll(allPoints, alli)
                {
                    allNearest[alli].proci = -1;
                    allNearest[alli].elementi = -1;
                    allNearest[alli].data = Foam::sqr(great);
                }
            }
            else
            {
                const labelList patchFaces(identity(pp.size()) + pp.start());

                treeBoundBox patchBb
                (
                    treeBoundBox(pp.points(), pp.meshPoints()).extend(1e-4)
                );

                indexedOctree<treeDataFace> boundaryTree
                (
                    treeDataFace    // all information needed to search faces
                    (
                        false,      // do not cache bb
                        mesh,
                        patchFaces  // boundary faces only
                    ),
                    patchBb,        // overall search domain
                    8,              // maxLevel
                    10,             // leafsize
                    3.0             // duplicity
                );

                forAll(allPoints, alli)
                {
                    const point& p = allPoints[alli];

                    const pointIndexHit pih =
                        boundaryTree.findNearest(p, magSqr(patchBb.span()));

                    if (pih.hit())
                    {
                        const point fc = pp[pih.index()].centre(pp.points());

                        allNearest[alli].proci = Pstream::myProcNo();
                        allNearest[alli].elementi = pih.index();
                        allNearest[alli].data = magSqr(fc - p);
                    }
                }
            }
            break;
        }

        case NEARESTFACE:
        {
            if (samplePatch().size() && samplePatch() != "none")
            {
                FatalErrorInFunction
                    << "No need to supply a patch name when in "
                    << sampleModeNames_[mode_] << " mode." << exit(FatalError);
            }

            //- Note: face-diagonal decomposition
            const meshSearchMeshObject& meshSearchEngine =
                meshSearchMeshObject::New(mesh);

            forAll(allPoints, alli)
            {
                const point& p = allPoints[alli];

                const label facei = meshSearchEngine.findNearestFace(p);

                if (facei != -1)
                {
                    const point& fc = mesh.faceCentres()[facei];

                    allNearest[alli].proci = Pstream::myProcNo();
                    allNearest[alli].elementi = facei;
                    allNearest[alli].data = magSqr(fc - p);
                }
            }
            break;
        }

        case NEARESTPATCHFACEAMI:
        {
            break;
        }
    }

    // Find nearest. Combine on master.
    Pstream::listCombineGather
    (
        allNearest,
        RemoteData<scalar>::smallestEqOp()
    );
    Pstream::listCombineScatter(allNearest);

    // Determine the number of missing samples (no reduction necessary as this
    // is computed from synchronised data)
    label nNotFound = 0;
    forAll(allPoints, alli)
    {
        if (allNearest[alli].proci == -1)
        {
            nNotFound ++;
        }
    }

    // If any points were not found within cells then re-search for them using
    // a nearest test, which should not fail. Warn that this is happening. If
    // any points were not found for some other method, then fail.
    if (nNotFound)
    {
        if (mode_ == NEARESTCELL)
        {
            WarningInFunction
                << "Did not find " << nNotFound
                << " out of " << returnReduce(patch_.size(), sumOp<label>())
                << " total samples. Sampling these on the nearest cell"
                << " centre instead." << endl << "On patch " << patch_.name()
                << " on region " << sampleRegion()
                << " with sample mode " << sampleModeNames_[mode_] << endl
                << "and offset mode " << offsetModeNames_[offsetMode_] << "."
                << endl;

            //- Note: face-diagonal decomposition
            const indexedOctree<Foam::treeDataCell>& tree = mesh.cellTree();

            forAll(allPoints, alli)
            {
                const point& p = allPoints[alli];

                if (allNearest[alli].proci == -1)
                {
                    const pointIndexHit pih = tree.findNearest(p, sqr(great));

                    allNearest[alli].proci = Pstream::myProcNo();
                    allNearest[alli].elementi = pih.index();
                    allNearest[alli].data = magSqr(pih.hitPoint() - p);
                }
            }

            Pstream::listCombineGather
            (
                allNearest,
                RemoteData<scalar>::smallestEqOp()
            );
            Pstream::listCombineScatter(allNearest);
        }
        else
        {
            FatalErrorInFunction
                << "Mapping failed for " << nl
                << "    patch: " << patch_.name() << nl
                << "    sampleRegion: " << sampleRegion() << nl
                << "    samplePatch: " << samplePatch() << nl
                << "    sampleMode: " << sampleModeNames_[mode_] << nl
                << "    offsetMode: " << offsetModeNames_[offsetMode_]
                << exit(FatalError);
        }
    }

    // Build lists of samples
    DynamicList<label> samplePatchGlobalFacesDyn;
    DynamicList<label> sampleIndicesDyn;
    forAll(allNearest, alli)
    {
        if (allNearest[alli].proci == Pstream::myProcNo())
        {
            samplePatchGlobalFacesDyn.append(alli);
            sampleIndicesDyn.append(allNearest[alli].elementi);
        }
    }
    sampleGlobalPatchFaces.transfer(samplePatchGlobalFacesDyn);
    sampleIndices.transfer(sampleIndicesDyn);
}


Foam::label Foam::mappedPatchBase::sampleSize() const
{
    switch (mode_)
    {
        case NEARESTCELL:
        {
            return sampleMesh().nCells();
        }
        case NEARESTPATCHFACE:
        {
            return samplePolyPatch().size();
        }
        case NEARESTPATCHFACEAMI:
        {
            return samplePolyPatch().size();
        }
        case NEARESTFACE:
        {
            return sampleMesh().nFaces() - sampleMesh().nInternalFaces();
        }
    }

    return -1;
}


void Foam::mappedPatchBase::calcMapping() const
{
    if (mapPtr_.valid())
    {
        FatalErrorInFunction
            << "Mapping already calculated" << exit(FatalError);
    }

    // Do a sanity check. Am I sampling my own patch? This only makes sense if
    // the position is transformed.
    if
    (
        mode_ == NEARESTPATCHFACE
     && sampleRegion() == patch_.boundaryMesh().mesh().name()
     && samplePatch() == patch_.name()
     && offsetMode_ == NONE
    )
    {
        FatalErrorInFunction
            << "Patch " << patch_.name() << " is sampling itself with no "
            << "offset. The patch face values are undefined."
            << exit(FatalError);

    }

    const globalIndex patchGlobalIndex(patch_.size());

    // Find processor and cell/face indices of samples
    labelList sampleGlobalPatchFaces, sampleIndices;
    findSamples(patchGlobalIndex, sampleGlobalPatchFaces, sampleIndices);

    // Construct distribution schedule
    List<Map<label>> compactMap;
    mapPtr_.reset
    (
        new distributionMap
        (
            patchGlobalIndex,
            sampleGlobalPatchFaces,
            compactMap
        )
    );
    const labelList oldToNew(move(sampleGlobalPatchFaces));
    const labelList oldSampleIndices(move(sampleIndices));

    // Construct input mapping for data to be distributed
    mapIndices_ = labelList(mapPtr_->constructSize(), -1);
    UIndirectList<label>(mapIndices_, oldToNew) = oldSampleIndices;

    // Reverse the map. This means the map is "forward" when going from samples
    // to the patch, which is logical.
    mapPtr_.reset
    (
        new distributionMap
        (
            patch_.size(),
            move(mapPtr_->constructMap()),
            move(mapPtr_->subMap())
        )
    );

    // Dump connecting lines for debugging
    if (debug && mode_ == NEARESTCELL)
    {
        OBJstream obj
        (
            patch_.name()
          + "_processor"
          + name(Pstream::myProcNo())
          + ".obj"
        );

        // Samples -> Patch
        pointField ccs(sampleMesh().cellCentres(), mapIndices_);
        mapPtr_->distribute(ccs);
        forAll(patch_, patchFacei)
        {
            const point& fc = patch_.faceCentres()[patchFacei];
            const point mid = 0.51*fc + 0.49*ccs[patchFacei];
            obj.write(linePointRef(fc, mid));
        }

        // Patch -> Samples
        pointField fcs(patch_.faceCentres());
        mapPtr_->reverseDistribute(mapIndices_.size(), fcs);
        forAll(mapIndices_, samplei)
        {
            const label celli = mapIndices_[samplei];
            if (celli == -1) continue;
            const point& cc = sampleMesh().cellCentres()[celli];
            const point mid = 0.51*cc + 0.49*fcs[samplei];
            obj.write(linePointRef(cc, mid));
        }
    }
}


void Foam::mappedPatchBase::calcAMI() const
{
    if (AMIPtr_.valid())
    {
        FatalErrorInFunction
            << "AMI already calculated" << exit(FatalError);
    }

    AMIPtr_.clear();

    // Get the projection surface, if any
    const word surfType(surfDict_.lookupOrDefault<word>("type", "none"));
    if (!surfPtr_.valid() && surfType != "none")
    {
        word surfName(surfDict_.lookupOrDefault("name", patch_.name()));

        const polyMesh& mesh = patch_.boundaryMesh().mesh();

        surfPtr_ =
            searchableSurface::New
            (
                surfType,
                IOobject
                (
                    surfName,
                    mesh.time().constant(),
                    searchableSurface::geometryDir(mesh.time()),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                surfDict_
            );
    }

    // Construct/apply AMI interpolation to determine addressing and weights
    AMIPtr_.reset
    (
        new AMIInterpolation
        (
            patch_,
            samplePolyPatch(), // nbrPatch0,
            surfPtr_,
            faceAreaIntersect::tmMesh,
            true,
            faceAreaWeightAMI::typeName,
            -1,
            AMIReverse_
        )
    );
}


// Hack to read old (List-based) format. See Field.C. The difference
// is only that in case of falling back to old format it expects a non-uniform
// list instead of a single vector.
Foam::tmp<Foam::pointField> Foam::mappedPatchBase::readListOrField
(
    const word& keyword,
    const dictionary& dict,
    const label size
)
{
    tmp<pointField> tfld(new pointField());
    pointField& fld = tfld.ref();

    if (size)
    {
        ITstream& is = dict.lookup(keyword);

        // Read first token
        token firstToken(is);

        if (firstToken.isWord())
        {
            if (firstToken.wordToken() == "uniform")
            {
                fld.setSize(size);
                fld = pTraits<vector>(is);
            }
            else if (firstToken.wordToken() == "nonuniform")
            {
                is >> static_cast<List<vector>&>(fld);
                if (fld.size() != size)
                {
                    FatalIOErrorInFunction
                    (
                        dict
                    )   << "size " << fld.size()
                        << " is not equal to the given value of " << size
                        << exit(FatalIOError);
                }
            }
            else
            {
                FatalIOErrorInFunction
                (
                    dict
                )   << "expected keyword 'uniform' or 'nonuniform', found "
                    << firstToken.wordToken()
                    << exit(FatalIOError);
            }
        }
        else
        {
            if (is.version() == 2.0)
            {
                IOWarningInFunction
                (
                    dict
                )   << "expected keyword 'uniform' or 'nonuniform', "
                       "assuming List format for backwards compatibility."
                       "Foam version 2.0." << endl;

                is.putBack(firstToken);
                is >> static_cast<List<vector>&>(fld);
            }
        }
    }
    return tfld;
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::mappedPatchBase::mappedPatchBase
(
    const polyPatch& pp
)
:
    patch_(pp),
    sampleRegion_(patch_.boundaryMesh().mesh().name()),
    sameRegion_(sampleRegion_ == patch_.boundaryMesh().mesh().name()),
    mode_(NEARESTPATCHFACE),
    samplePatch_(word::null),
    coupleGroup_(),
    offsetMode_(NONE),
    distance_(NaN),
    offset_(vector::uniform(NaN)),
    mapPtr_(nullptr),
    mapIndices_(),
    AMIPtr_(nullptr),
    AMIReverse_(false),
    surfPtr_(nullptr),
    surfDict_(fileName("surface"))
{}


Foam::mappedPatchBase::mappedPatchBase
(
    const polyPatch& pp,
    const word& sampleRegion,
    const sampleMode mode,
    const word& samplePatch
)
:
    patch_(pp),
    sampleRegion_(sampleRegion),
    sameRegion_(sampleRegion_ == patch_.boundaryMesh().mesh().name()),
    mode_(mode),
    samplePatch_(samplePatch),
    coupleGroup_(),
    offsetMode_(NONE),
    distance_(NaN),
    offset_(vector::uniform(NaN)),
    mapPtr_(nullptr),
    mapIndices_(),
    AMIPtr_(nullptr),
    AMIReverse_(false),
    surfPtr_(nullptr),
    surfDict_(fileName("surface"))
{}


Foam::mappedPatchBase::mappedPatchBase
(
    const polyPatch& pp,
    const dictionary& dict
)
:
    patch_(pp),
    sampleRegion_(dict.lookupOrDefault<word>("sampleRegion", word::null)),
    sameRegion_(sampleRegion_ == patch_.boundaryMesh().mesh().name()),
    mode_(sampleModeNames_.read(dict.lookup("sampleMode"))),
    samplePatch_(dict.lookupOrDefault<word>("samplePatch", word::null)),
    coupleGroup_(dict),
    offsetMode_(readOffsetMode(dict)),
    distance_
    (
        offsetMode_ == NORMAL
      ? dict.lookup<scalar>("distance")
      : NaN
    ),
    offset_
    (
        offsetMode_ == DIRECTION
      ? dict.lookup<vector>("offset")
      : vector::uniform(NaN)
    ),
    mapPtr_(nullptr),
    mapIndices_(),
    AMIPtr_(nullptr),
    AMIReverse_(dict.lookupOrDefault<bool>("flipNormals", false)),
    surfPtr_(nullptr),
    surfDict_(dict.subOrEmptyDict("surface"))
{
    if (!coupleGroup_.valid() && sampleRegion_.empty())
    {
        // If no coupleGroup and no sampleRegion then assume the local region
        sampleRegion_ = patch_.boundaryMesh().mesh().name();
        sameRegion_ = true;
    }
}


Foam::mappedPatchBase::mappedPatchBase
(
    const polyPatch& pp,
    const mappedPatchBase& mpb
)
:
    patch_(pp),
    sampleRegion_(mpb.sampleRegion_),
    sameRegion_(mpb.sameRegion_),
    mode_(mpb.mode_),
    samplePatch_(mpb.samplePatch_),
    coupleGroup_(mpb.coupleGroup_),
    offsetMode_(mpb.offsetMode_),
    distance_(mpb.distance_),
    offset_(mpb.offset_),
    mapPtr_(nullptr),
    mapIndices_(),
    AMIPtr_(nullptr),
    AMIReverse_(mpb.AMIReverse_),
    surfPtr_(nullptr),
    surfDict_(mpb.surfDict_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mappedPatchBase::~mappedPatchBase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::polyMesh& Foam::mappedPatchBase::sampleMesh() const
{
    return patch_.boundaryMesh().mesh().time().lookupObject<polyMesh>
    (
        sampleRegion()
    );
}


const Foam::polyPatch& Foam::mappedPatchBase::samplePolyPatch() const
{
    const polyMesh& nbrMesh = sampleMesh();

    const label patchi = nbrMesh.boundaryMesh().findPatchID(samplePatch());

    if (patchi == -1)
    {
        FatalErrorInFunction
            << "Cannot find patch " << samplePatch()
            << " in region " << sampleRegion_ << endl
            << "Valid patches are " << nbrMesh.boundaryMesh().names()
            << exit(FatalError);
    }

    return nbrMesh.boundaryMesh()[patchi];
}


Foam::tmp<Foam::pointField> Foam::mappedPatchBase::samplePoints() const
{
    const polyMesh& mesh = patch_.boundaryMesh().mesh();

    // Force construction of min-tet decomp
    (void)mesh.tetBasePtIs();

    // Allocate the result
    tmp<pointField> tresult(new pointField(patch_.size()));
    pointField& result = tresult.ref();

    // Compute the face points. Assume the face centre to begin with, and then
    // try and find a better match on complex faces by doing a ray intersection
    // with a triangulation of the face. This triangulation matches the
    // tetrahedralisation used for cell searches. This maximises the chances
    // that the point will be successfully found in a cell.
    forAll(patch_, patchFacei)
    {
        const pointField& ps = mesh.points();

        const face& f = patch_[patchFacei];
        const point& fc = patch_.faceCentres()[patchFacei];

        result[patchFacei] = fc;

        if (f.size() > 3)
        {
            const label facei = patch_.start() + patchFacei;
            const label celli = patch_.faceCells()[patchFacei];

            const point& cc = mesh.cellCentres()[celli];
            const vector d = fc - cc;

            const label faceBasePtI = mesh.tetBasePtIs()[facei];

            for (label tetPtI = 1; tetPtI < f.size() - 1; tetPtI ++)
            {
                const label facePtI = (tetPtI + faceBasePtI) % f.size();
                const label otherFacePtI = f.fcIndex(facePtI);

                const triPointRef tri
                (
                    ps[f[faceBasePtI]],
                    ps[f[facePtI]],
                    ps[f[otherFacePtI]]
                );

                const pointHit hitInfo =
                    tri.intersection(cc, d, intersection::algorithm::halfRay);

                if (hitInfo.hit() && hitInfo.distance() > 0)
                {
                    result[patchFacei] = hitInfo.hitPoint();
                    break;
                }
            }
        }
    }

    // Apply offset to get sample points
    switch (offsetMode_)
    {
        case NONE:
            break;
        case NORMAL:
            result += distance_*patch_.faceNormals();
            break;
        case DIRECTION:
            result += offset_;
            break;
    }

    return tresult;
}


void Foam::mappedPatchBase::clearOut()
{
    mapPtr_.clear();
    mapIndices_.clear();
    AMIPtr_.clear();
    surfPtr_.clear();
}


void Foam::mappedPatchBase::write(Ostream& os) const
{
    writeEntry(os, "sampleMode", sampleModeNames_[mode_]);
    if (!sampleRegion_.empty())
    {
        writeEntry(os, "sampleRegion", sampleRegion_);
    }
    if (!samplePatch_.empty())
    {
        writeEntry(os, "samplePatch", samplePatch_);
    }
    coupleGroup_.write(os);

    writeEntry(os, "offsetMode", offsetModeNames_[offsetMode_]);

    switch (offsetMode_)
    {
        case NONE:
            break;
        case NORMAL:
            writeEntry(os, "distance", distance_);
            break;
        case DIRECTION:
            writeEntry(os, "offset", offset_);
            break;
    }

    if (mode_ == NEARESTPATCHFACEAMI)
    {
        if (AMIReverse_)
        {
            writeEntry(os, "flipNormals", AMIReverse_);
        }

        if (!surfDict_.empty())
        {
            writeKeyword(os, surfDict_.dictName());
            os  << surfDict_;
        }
    }
}


// ************************************************************************* //
