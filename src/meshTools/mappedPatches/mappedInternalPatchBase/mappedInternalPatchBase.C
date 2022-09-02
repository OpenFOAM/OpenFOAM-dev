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

#include "mappedInternalPatchBase.H"
#include "SubField.H"
#include "Time.H"
#include "triPointRef.H"
#include "treeDataCell.H"
#include "indexedOctree.H"
#include "globalIndex.H"
#include "RemoteData.H"
#include "OBJstream.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mappedInternalPatchBase, 0);

    template<>
    const char* Foam::NamedEnum
    <
        Foam::mappedInternalPatchBase::offsetMode,
        2
    >::names[] =
    {
        "normal",
        "direction"
    };
}


const Foam::NamedEnum<Foam::mappedInternalPatchBase::offsetMode, 2>
    Foam::mappedInternalPatchBase::offsetModeNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::mappedInternalPatchBase::offsetMode
Foam::mappedInternalPatchBase::readOffsetMode
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

        if (haveDistance == haveOffset)
        {
            // Error. Demand "offsetMode" setting to disambiguate.
            return offsetModeNames_.read(dict.lookup("offsetMode"));
        }
        else if (haveDistance)
        {
            return NORMAL;
        }
        else
        {
            return DIRECTION;
        }
    }
}


void Foam::mappedInternalPatchBase::calcMapping() const
{
    if (mapPtr_.valid())
    {
        FatalErrorInFunction
            << "Mapping already calculated" << exit(FatalError);
    }

    const polyMesh& nbrMesh = this->nbrMesh();

    const globalIndex patchGlobalIndex(patch_.size());

    // Find processor and cell/face indices of samples
    labelList sampleGlobalPatchFaces, sampleIndices;
    {
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

        // List of possibly remote sampling cells
        List<RemoteData<scalar>> allNearest(patchGlobalIndex.size());

        // Find containing cell for every sampling point
        const indexedOctree<Foam::treeDataCell>& tree = nbrMesh.cellTree();
        forAll(allPoints, alli)
        {
            const point& p = allPoints[alli];

            const label celli = tree.findInside(p);

            if (celli != -1)
            {
                const point& cc = nbrMesh.cellCentres()[celli];

                allNearest[alli].proci = Pstream::myProcNo();
                allNearest[alli].elementi = celli;
                allNearest[alli].data = magSqr(cc - p);
            }
        }

        // Find nearest. Combine on master.
        Pstream::listCombineGather
        (
            allNearest,
            RemoteData<scalar>::smallestEqOp()
        );
        Pstream::listCombineScatter(allNearest);

        // Determine the number of missing samples (no reduction necessary as
        // this is computed from synchronised data)
        label nNotFound = 0;
        forAll(allPoints, alli)
        {
            if (allNearest[alli].proci == -1)
            {
                nNotFound ++;
            }
        }

        // If any points were not found within cells then re-search for them
        // using a nearest test, which should not fail. Warn that this is
        // happening. If any points were not found for some other method, then
        // fail.
        if (nNotFound)
        {
            WarningInFunction
                << "Did not find a containing cell for " << nNotFound
                << " out of " << returnReduce(patch_.size(), sumOp<label>())
                << " total faces. Using nearest cell for these faces instead."
                << endl << "On patch " << patch_.name() << " on region "
                << nbrRegionName() << " with offset mode "
                << offsetModeNames_[offsetMode_] << "." << endl;

            const indexedOctree<Foam::treeDataCell>& tree = nbrMesh.cellTree();
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
    cellIndices_ = labelList(mapPtr_->constructSize(), -1);
    UIndirectList<label>(cellIndices_, oldToNew) = oldSampleIndices;

    // Reverse the map. This means the map is "forward" when going from cells
    // to this patch, which is logical.
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
    if (debug)
    {
        OBJstream obj
        (
            patch_.name()
          + "_processor"
          + name(Pstream::myProcNo())
          + ".obj"
        );

        // Cells -> Patch
        pointField ccs(nbrMesh.cellCentres(), cellIndices_);
        mapPtr_->distribute(ccs);
        forAll(patch_, patchFacei)
        {
            const point& fc = patch_.faceCentres()[patchFacei];
            const point mid = 0.51*fc + 0.49*ccs[patchFacei];
            obj.write(linePointRef(fc, mid));
        }

        // Patch -> Cells
        pointField fcs(patch_.faceCentres());
        mapPtr_->reverseDistribute(cellIndices_.size(), fcs);
        forAll(cellIndices_, i)
        {
            const label celli = cellIndices_[i];
            if (celli == -1) continue;
            const point& cc = nbrMesh.cellCentres()[celli];
            const point mid = 0.51*cc + 0.49*fcs[i];
            obj.write(linePointRef(cc, mid));
        }
    }
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::mappedInternalPatchBase::mappedInternalPatchBase
(
    const polyPatch& pp
)
:
    patch_(pp),
    nbrRegionName_(patch_.boundaryMesh().mesh().name()),
    offsetMode_(NORMAL),
    distance_(NaN),
    offset_(vector::uniform(NaN)),
    mapPtr_(nullptr),
    cellIndices_()
{}


Foam::mappedInternalPatchBase::mappedInternalPatchBase
(
    const polyPatch& pp,
    const word& neighbourRegion
)
:
    patch_(pp),
    nbrRegionName_(neighbourRegion),
    offsetMode_(NORMAL),
    distance_(NaN),
    offset_(vector::uniform(NaN)),
    mapPtr_(nullptr),
    cellIndices_()
{}


Foam::mappedInternalPatchBase::mappedInternalPatchBase
(
    const polyPatch& pp,
    const dictionary& dict
)
:
    patch_(pp),
    nbrRegionName_
    (
        dict.lookupOrDefaultBackwardsCompatible<word>
        (
            {"neighbourRegion", "sampleRegion"},
            patch_.boundaryMesh().mesh().name()
        )
    ),
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
    cellIndices_()
{}


Foam::mappedInternalPatchBase::mappedInternalPatchBase
(
    const polyPatch& pp,
    const mappedInternalPatchBase& mpb
)
:
    patch_(pp),
    nbrRegionName_(mpb.nbrRegionName_),
    offsetMode_(mpb.offsetMode_),
    distance_(mpb.distance_),
    offset_(mpb.offset_),
    mapPtr_(nullptr),
    cellIndices_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mappedInternalPatchBase::~mappedInternalPatchBase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::polyMesh& Foam::mappedInternalPatchBase::nbrMesh() const
{
    return patch_.boundaryMesh().mesh().time().lookupObject<polyMesh>
    (
        nbrRegionName()
    );
}


Foam::tmp<Foam::pointField> Foam::mappedInternalPatchBase::samplePoints() const
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
        case NORMAL:
            result += distance_*patch_.faceNormals();
            break;
        case DIRECTION:
            result += offset_;
            break;
    }

    return tresult;
}


void Foam::mappedInternalPatchBase::clearOut()
{
    mapPtr_.clear();
    cellIndices_.clear();
}


bool Foam::mappedInternalPatchBase::specified(const dictionary& dict)
{
    return
        dict.found("neighbourRegion")
     || dict.found("offsetMode")
     || dict.found("distance")
     || dict.found("offset");
}


void Foam::mappedInternalPatchBase::write(Ostream& os) const
{
    writeEntryIfDifferent
    (
        os,
        "neighbourRegion",
        patch_.boundaryMesh().mesh().name(),
        nbrRegionName_
    );

    writeEntry(os, "offsetMode", offsetModeNames_[offsetMode_]);

    switch (offsetMode_)
    {
        case NORMAL:
            writeEntry(os, "distance", distance_);
            break;
        case DIRECTION:
            writeEntry(os, "offset", offset_);
            break;
    }
}


// ************************************************************************* //
