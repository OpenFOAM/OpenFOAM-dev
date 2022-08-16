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
#include "addToRunTimeSelectionTable.H"
#include "ListListOps.H"
#include "meshSearchMeshObject.H"
#include "meshTools.H"
#include "OFstream.H"
#include "Random.H"
#include "treeDataFace.H"
#include "treeDataPoint.H"
#include "indexedOctree.H"
#include "polyMesh.H"
#include "polyPatch.H"
#include "Time.H"
#include "distributionMap.H"
#include "SubField.H"
#include "triPointRef.H"
#include "syncTools.H"
#include "treeDataCell.H"
#include "DynamicField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mappedPatchBase, 0);

    template<>
    const char* Foam::NamedEnum
    <
        Foam::mappedPatchBase::sampleMode,
        5
    >::names[] =
    {
        "nearestCell",
        "nearestPatchFace",
        "nearestPatchFaceAMI",
        "nearestFace",
        "nearestOnlyCell"
    };

    template<>
    const char* Foam::NamedEnum
    <
        Foam::mappedPatchBase::offsetMode,
        3
    >::names[] =
    {
        "uniform",
        "nonuniform",
        "normal"
    };
}


const Foam::NamedEnum<Foam::mappedPatchBase::sampleMode, 5>
    Foam::mappedPatchBase::sampleModeNames_;

const Foam::NamedEnum<Foam::mappedPatchBase::offsetMode, 3>
    Foam::mappedPatchBase::offsetModeNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::mappedPatchBase::collectSamples
(
    pointField& samples,
    labelList& patchFaceProcs,
    labelList& patchFaces
) const
{
    // Collect all sample points and the faces they come from.

    {
        List<pointField> globalSamples(Pstream::nProcs());
        globalSamples[Pstream::myProcNo()] = samplePoints();
        Pstream::gatherList(globalSamples);
        Pstream::scatterList(globalSamples);

        samples = ListListOps::combine<pointField>
        (
            globalSamples,
            accessOp<pointField>()
        );
    }

    {
        labelListList globalFaces(Pstream::nProcs());
        globalFaces[Pstream::myProcNo()] = identity(patch_.size());
        Pstream::gatherList(globalFaces);
        Pstream::scatterList(globalFaces);

        patchFaces = ListListOps::combine<labelList>
        (
            globalFaces,
            accessOp<labelList>()
        );
    }

    {
        labelList nPerProc(Pstream::nProcs());
        nPerProc[Pstream::myProcNo()] = patch_.size();
        Pstream::gatherList(nPerProc);
        Pstream::scatterList(nPerProc);

        patchFaceProcs.setSize(patchFaces.size());

        label sampleI = 0;
        forAll(nPerProc, proci)
        {
            for (label i = 0; i < nPerProc[proci]; i++)
            {
                patchFaceProcs[sampleI++] = proci;
            }
        }
    }
}


void Foam::mappedPatchBase::findSamples
(
    const sampleMode mode,
    const pointField& samples,
    labelList& sampleProcs,
    labelList& sampleIndices
) const
{
    // Lookup the correct region
    const polyMesh& mesh = sampleMesh();

    // All the info for nearest. Construct to miss
    List<nearInfo> nearest(samples.size());

    switch (mode)
    {
        case NEARESTCELL:
        {
            if (samplePatch_.size() && samplePatch_ != "none")
            {
                FatalErrorInFunction
                    << "No need to supply a patch name when in "
                    << sampleModeNames_[mode] << " mode." << exit(FatalError);
            }

            //- Note: face-diagonal decomposition
            const indexedOctree<Foam::treeDataCell>& tree = mesh.cellTree();

            forAll(samples, sampleI)
            {
                const point& sample = samples[sampleI];

                label celli = tree.findInside(sample);

                if (celli == -1)
                {
                    nearest[sampleI].second().first() = Foam::sqr(great);
                    nearest[sampleI].second().second() = Pstream::myProcNo();
                }
                else
                {
                    const point& cc = mesh.cellCentres()[celli];

                    nearest[sampleI].first() = pointIndexHit
                    (
                        true,
                        cc,
                        celli
                    );
                    nearest[sampleI].second().first() = magSqr(cc-sample);
                    nearest[sampleI].second().second() = Pstream::myProcNo();
                }
            }
            break;
        }

        case NEARESTONLYCELL:
        {
            if (samplePatch_.size() && samplePatch_ != "none")
            {
                FatalErrorInFunction
                    << "No need to supply a patch name when in "
                    << sampleModeNames_[mode] << " mode." << exit(FatalError);
            }

            //- Note: face-diagonal decomposition
            const indexedOctree<Foam::treeDataCell>& tree = mesh.cellTree();

            forAll(samples, sampleI)
            {
                const point& sample = samples[sampleI];

                nearest[sampleI].first() = tree.findNearest(sample, sqr(great));
                nearest[sampleI].second().first() = magSqr
                (
                    nearest[sampleI].first().hitPoint()
                   -sample
                );
                nearest[sampleI].second().second() = Pstream::myProcNo();
            }
            break;
        }

        case NEARESTPATCHFACE:
        {
            const polyPatch& pp = samplePolyPatch();

            if (pp.empty())
            {
                forAll(samples, sampleI)
                {
                    nearest[sampleI].second().first() = Foam::sqr(great);
                    nearest[sampleI].second().second() = Pstream::myProcNo();
                }
            }
            else
            {
                // patch faces
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

                forAll(samples, sampleI)
                {
                    const point& sample = samples[sampleI];

                    pointIndexHit& nearInfo = nearest[sampleI].first();
                    nearInfo = boundaryTree.findNearest
                    (
                        sample,
                        magSqr(patchBb.span())
                    );

                    if (!nearInfo.hit())
                    {
                        nearest[sampleI].second().first() = Foam::sqr(great);
                        nearest[sampleI].second().second() =
                            Pstream::myProcNo();
                    }
                    else
                    {
                        point fc(pp[nearInfo.index()].centre(pp.points()));
                        nearInfo.setPoint(fc);
                        nearest[sampleI].second().first() = magSqr(fc-sample);
                        nearest[sampleI].second().second() =
                            Pstream::myProcNo();
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
                    << sampleModeNames_[mode] << " mode." << exit(FatalError);
            }

            //- Note: face-diagonal decomposition
            const meshSearchMeshObject& meshSearchEngine =
                meshSearchMeshObject::New(mesh);

            forAll(samples, sampleI)
            {
                const point& sample = samples[sampleI];

                label facei = meshSearchEngine.findNearestFace(sample);

                if (facei == -1)
                {
                    nearest[sampleI].second().first() = Foam::sqr(great);
                    nearest[sampleI].second().second() = Pstream::myProcNo();
                }
                else
                {
                    const point& fc = mesh.faceCentres()[facei];

                    nearest[sampleI].first() = pointIndexHit
                    (
                        true,
                        fc,
                        facei
                    );
                    nearest[sampleI].second().first() = magSqr(fc-sample);
                    nearest[sampleI].second().second() = Pstream::myProcNo();
                }
            }
            break;
        }

        case NEARESTPATCHFACEAMI:
        {
            // nothing to do here
            return;
        }

        default:
        {
            FatalErrorInFunction
                << "problem." << abort(FatalError);
        }
    }

    // Find nearest. Combine on master.
    Pstream::listCombineGather(nearest, nearestEqOp());
    Pstream::listCombineScatter(nearest);

    if (debug)
    {
        InfoInFunction
            << "mesh " << sampleRegion() << " : " << endl;

        forAll(nearest, sampleI)
        {
            label proci = nearest[sampleI].second().second();
            label localI = nearest[sampleI].first().index();

            Info<< "    " << sampleI << " coord:"<< samples[sampleI]
                << " found on processor:" << proci
                << " in local cell/face/point:" << localI
                << " with location:" << nearest[sampleI].first().rawPoint()
                << endl;
        }
    }

    // Convert back into proc+local index
    sampleProcs.setSize(samples.size());
    sampleIndices.setSize(samples.size());

    forAll(nearest, sampleI)
    {
        if (!nearest[sampleI].first().hit())
        {
            sampleProcs[sampleI] = -1;
            sampleIndices[sampleI] = -1;
        }
        else
        {
            sampleProcs[sampleI] = nearest[sampleI].second().second();
            sampleIndices[sampleI] = nearest[sampleI].first().index();
        }
    }
}


Foam::label Foam::mappedPatchBase::sampleSize() const
{
    switch (mode_)
    {
        case NEARESTPATCHFACEAMI:
        {
            return samplePolyPatch().size();
        }
        case NEARESTCELL:
        {
            return sampleMesh().nCells();
        }
        case NEARESTPATCHFACE:
        {
            return samplePolyPatch().size();
        }
        case NEARESTFACE:
        {
            return sampleMesh().nFaces() - sampleMesh().nInternalFaces();
        }
        default:
        {
            FatalErrorInFunction
                << "problem." << abort(FatalError);
            return -1;
        }
    }
}


void Foam::mappedPatchBase::calcMapping() const
{
    if (mapPtr_.valid())
    {
        FatalErrorInFunction
            << "Mapping already calculated" << exit(FatalError);
    }

    /*
    // Do a sanity check. Am I sampling my own patch? This only makes sense if
    // the position is transformed.
    if
    (
        mode_ == NEARESTPATCHFACE
     && sampleRegion() == patch_.boundaryMesh().mesh().name()
     && samplePatch() == patch_.name()
     && !transform_.transformsPosition()
    )
    {
        FatalErrorInFunction
            << "Patch " << patch_.name() << " is sampling itself with no "
            << "transformation. The patch face values are undefined."
            << exit(FatalError);

    }
    */

    // Get global list of all samples and the processor and face they come from.
    pointField samples;
    labelList patchFaceProcs;
    labelList patchFaces;
    collectSamples(samples, patchFaceProcs, patchFaces);

    // Find processor and cell/face samples are in and actual location.
    labelList sampleProcs;
    labelList sampleIndices;
    findSamples(mode_, samples, sampleProcs, sampleIndices);

    // Check for samples that were not found. This will only happen for
    // NEARESTCELL since this finds a cell containing a location.
    if (mode_ == NEARESTCELL)
    {
        const label nNotFound =
            returnReduce(count(sampleProcs, -1), sumOp<label>());

        if (nNotFound > 0)
        {
            WarningInFunction
                << "Did not find " << nNotFound
                << " out of " << sampleProcs.size() << " total samples."
                << " Sampling these on the nearest cell centre instead." << endl
                << "On patch " << patch_.name()
                << " on region " << sampleRegion()
                << " with sample mode " << sampleModeNames_[mode_] << endl
                << "and offset mode " << offsetModeNames_[offsetMode_] << "."
                << endl;

            // Collect the samples that cannot be found
            DynamicList<label> subMap;
            DynamicField<point> subSamples;
            forAll(sampleProcs, sampleI)
            {
                if (sampleProcs[sampleI] == -1)
                {
                    subMap.append(sampleI);
                    subSamples.append(samples[sampleI]);
                }
            }

            // And re-search for pure nearest (should not fail)
            labelList subSampleProcs;
            labelList subSampleIndices;
            findSamples
            (
                NEARESTONLYCELL,
                subSamples,
                subSampleProcs,
                subSampleIndices
            );

            // Insert
            UIndirectList<label>(sampleProcs, subMap) = subSampleProcs;
            UIndirectList<label>(sampleIndices, subMap) = subSampleIndices;
        }
    }

    // Check for mapping failure
    if (findIndex(sampleProcs, -1) != -1)
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

    // Determine schedule.
    mapPtr_.reset(new distributionMap(sampleProcs, patchFaceProcs));

    // Rework the schedule from indices into samples to cell data to send,
    // face data to receive.
    labelListList& subMap = mapPtr_().subMap();
    labelListList& constructMap = mapPtr_().constructMap();
    forAll(subMap, proci)
    {
        subMap[proci] = UIndirectList<label>
        (
            sampleIndices,
            subMap[proci]
        );
        constructMap[proci] = UIndirectList<label>
        (
            patchFaces,
            constructMap[proci]
        );
    }

    // Redo constructSize
    mapPtr_().constructSize() = patch_.size();

    if (debug)
    {
        // Check that all elements get a value.
        PackedBoolList used(patch_.size(), false);

        forAll(constructMap, proci)
        {
            const labelList& map = constructMap[proci];

            forAll(map, i)
            {
                label facei = map[i];

                if (!used[facei])
                {
                    used[facei] = true;
                }
                else
                {
                    FatalErrorInFunction
                        << "On patch " << patch_.name()
                        << " patchface " << facei
                        << " is assigned to more than once."
                        << abort(FatalError);
                }
            }
        }

        forAll(used, facei)
        {
            if (used[facei] == 0)
            {
                FatalErrorInFunction
                    << "On patch " << patch_.name()
                    << " patchface " << facei
                    << " is never assigned to."
                    << abort(FatalError);
            }
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

    if (debug)
    {
        const polyPatch& nbr = samplePolyPatch();

        pointField nbrPoints(nbr.localPoints());

        OFstream os(patch_.name() + "_neighbourPatch-org.obj");
        meshTools::writeOBJ(os, samplePolyPatch().localFaces(), nbrPoints);

        // transform neighbour patch to local system
        primitivePatch nbrPatch0
        (
            SubList<face>
            (
                nbr.localFaces(),
                nbr.size()
            ),
            nbrPoints
        );

        OFstream osN(patch_.name() + "_neighbourPatch-trans.obj");
        meshTools::writeOBJ(osN, nbrPatch0, nbrPoints);

        OFstream osO(patch_.name() + "_ownerPatch.obj");
        meshTools::writeOBJ(osO, patch_.localFaces(), patch_.localPoints());
    }

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
    mode_(NEARESTPATCHFACE),
    samplePatch_(""),
    coupleGroup_(),
    offsetMode_(UNIFORM),
    offset_(Zero),
    offsets_(pp.size(), offset_),
    distance_(0),
    sameRegion_(sampleRegion_ == patch_.boundaryMesh().mesh().name()),
    mapPtr_(nullptr),
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
    const word& samplePatch,
    const vectorField& offsets
)
:
    patch_(pp),
    sampleRegion_(sampleRegion),
    mode_(mode),
    samplePatch_(samplePatch),
    coupleGroup_(),
    offsetMode_(NONUNIFORM),
    offset_(Zero),
    offsets_(offsets),
    distance_(0),
    sameRegion_(sampleRegion_ == patch_.boundaryMesh().mesh().name()),
    mapPtr_(nullptr),
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
    const word& samplePatch,
    const vector& offset
)
:
    patch_(pp),
    sampleRegion_(sampleRegion),
    mode_(mode),
    samplePatch_(samplePatch),
    coupleGroup_(),
    offsetMode_(UNIFORM),
    offset_(offset),
    offsets_(0),
    distance_(0),
    sameRegion_(sampleRegion_ == patch_.boundaryMesh().mesh().name()),
    mapPtr_(nullptr),
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
    const word& samplePatch,
    const scalar distance
)
:
    patch_(pp),
    sampleRegion_(sampleRegion),
    mode_(mode),
    samplePatch_(samplePatch),
    coupleGroup_(),
    offsetMode_(NORMAL),
    offset_(Zero),
    offsets_(0),
    distance_(distance),
    sameRegion_(sampleRegion_ == patch_.boundaryMesh().mesh().name()),
    mapPtr_(nullptr),
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
    sampleRegion_(dict.lookupOrDefault<word>("sampleRegion", "")),
    mode_(sampleModeNames_.read(dict.lookup("sampleMode"))),
    samplePatch_(dict.lookupOrDefault<word>("samplePatch", "")),
    coupleGroup_(dict),
    offsetMode_(UNIFORM),
    offset_(Zero),
    offsets_(0),
    distance_(0.0),
    sameRegion_(sampleRegion_ == patch_.boundaryMesh().mesh().name()),
    mapPtr_(nullptr),
    AMIPtr_(nullptr),
    AMIReverse_(dict.lookupOrDefault<bool>("flipNormals", false)),
    surfPtr_(nullptr),
    surfDict_(dict.subOrEmptyDict("surface"))
{
    if (!coupleGroup_.valid())
    {
        if (sampleRegion_.empty())
        {
            // If no coupleGroup and no sampleRegion assume local region
            sampleRegion_ = patch_.boundaryMesh().mesh().name();
            sameRegion_ = true;
        }
    }

    if (dict.found("offsetMode"))
    {
        offsetMode_ = offsetModeNames_.read(dict.lookup("offsetMode"));

        switch(offsetMode_)
        {
            case UNIFORM:
            {
                offset_ = point(dict.lookup("offset"));
            }
            break;

            case NONUNIFORM:
            {
                offsets_ = readListOrField("offsets", dict, patch_.size());
            }
            break;

            case NORMAL:
            {
                distance_ = dict.lookup<scalar>("distance");
            }
            break;
        }
    }
    else if (dict.found("offset"))
    {
        offsetMode_ = UNIFORM;
        offset_ = point(dict.lookup("offset"));
    }
    else if (dict.found("offsets"))
    {
        offsetMode_ = NONUNIFORM;
        offsets_ = readListOrField("offsets", dict, patch_.size());
    }
    else if (mode_ != NEARESTPATCHFACE && mode_ != NEARESTPATCHFACEAMI)
    {
        FatalIOErrorInFunction(dict)
            << "Please supply the offsetMode as one of "
            << NamedEnum<offsetMode, 3>::words()
            << exit(FatalIOError);
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
    mode_(mpb.mode_),
    samplePatch_(mpb.samplePatch_),
    coupleGroup_(mpb.coupleGroup_),
    offsetMode_(mpb.offsetMode_),
    offset_(mpb.offset_),
    offsets_(mpb.offsets_),
    distance_(mpb.distance_),
    sameRegion_(mpb.sameRegion_),
    mapPtr_(nullptr),
    AMIPtr_(nullptr),
    AMIReverse_(mpb.AMIReverse_),
    surfPtr_(nullptr),
    surfDict_(mpb.surfDict_)
{}


Foam::mappedPatchBase::mappedPatchBase
(
    const polyPatch& pp,
    const mappedPatchBase& mpb,
    const labelUList& mapAddressing
)
:
    patch_(pp),
    sampleRegion_(mpb.sampleRegion_),
    mode_(mpb.mode_),
    samplePatch_(mpb.samplePatch_),
    coupleGroup_(mpb.coupleGroup_),
    offsetMode_(mpb.offsetMode_),
    offset_(mpb.offset_),
    offsets_
    (
        offsetMode_ == NONUNIFORM
      ? vectorField(mpb.offsets_, mapAddressing)
      : vectorField(0)
    ),
    distance_(mpb.distance_),
    sameRegion_(mpb.sameRegion_),
    mapPtr_(nullptr),
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
        case UNIFORM:
        {
            result += offset_;
            break;
        }
        case NONUNIFORM:
        {
            result += offsets_;
            break;
        }
        case NORMAL:
        {
            result += distance_*patch_.faceNormals();
            break;
        }
    }

    return tresult;
}


void Foam::mappedPatchBase::clearOut()
{
    mapPtr_.clear();
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

    if
    (
        offsetMode_ == UNIFORM
     && offset_ == vector::zero
     && (mode_ == NEARESTPATCHFACE || mode_ == NEARESTPATCHFACEAMI)
    )
    {
        // Collocated mode. No need to write offset data
    }
    else
    {
        writeEntry(os, "offsetMode", offsetModeNames_[offsetMode_]);

        switch (offsetMode_)
        {
            case UNIFORM:
            {
                writeEntry(os, "offset", offset_);
                break;
            }
            case NONUNIFORM:
            {
                writeEntry(os, "offsets", offsets_);
                break;
            }
            case NORMAL:
            {
                writeEntry(os, "distance", distance_);
                break;
            }
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
}


// ************************************************************************* //
