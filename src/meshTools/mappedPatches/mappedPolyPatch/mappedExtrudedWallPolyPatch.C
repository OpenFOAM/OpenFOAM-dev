/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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

#include "mappedExtrudedWallPolyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "LayerInfoData.H"
#include "FaceCellWave.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mappedExtrudedWallPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, mappedExtrudedWallPolyPatch, word);
    addToRunTimeSelectionTable
    (
        polyPatch,
        mappedExtrudedWallPolyPatch,
        dictionary
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::mappedExtrudedWallPolyPatch::initCalcGeometry(PstreamBuffers& pBufs)
{
    mappedWallPolyPatch::initCalcGeometry(pBufs);
}


void Foam::mappedExtrudedWallPolyPatch::calcGeometry(PstreamBuffers& pBufs)
{
    mappedWallPolyPatch::calcGeometry(pBufs);
    samplePointsPtr_.clear();
}


void Foam::mappedExtrudedWallPolyPatch::initMovePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    mappedWallPolyPatch::initMovePoints(pBufs, p);
}


void Foam::mappedExtrudedWallPolyPatch::movePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    mappedWallPolyPatch::movePoints(pBufs, p);
    samplePointsPtr_.clear();
}


void Foam::mappedExtrudedWallPolyPatch::initTopoChange(PstreamBuffers& pBufs)
{
    mappedWallPolyPatch::initTopoChange(pBufs);
}


void Foam::mappedExtrudedWallPolyPatch::topoChange(PstreamBuffers& pBufs)
{
    mappedWallPolyPatch::topoChange(pBufs);
    samplePointsPtr_.clear();
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::mappedExtrudedWallPolyPatch::mappedExtrudedWallPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    mappedWallPolyPatch(name, size, start, index, bm, patchType),
    bottomPatch_(word::null)
{}


Foam::mappedExtrudedWallPolyPatch::mappedExtrudedWallPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    mappedWallPolyPatch(name, dict, index, bm, patchType),
    bottomPatch_(dict.lookupOrDefault<word>("bottomPatch", word::null))
{}


Foam::mappedExtrudedWallPolyPatch::mappedExtrudedWallPolyPatch
(
    const mappedExtrudedWallPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    mappedWallPolyPatch(pp, bm),
    bottomPatch_(pp.bottomPatch_)
{}


Foam::mappedExtrudedWallPolyPatch::mappedExtrudedWallPolyPatch
(
    const mappedExtrudedWallPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    mappedWallPolyPatch(pp, bm, index, newSize, newStart),
    bottomPatch_(pp.bottomPatch_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mappedExtrudedWallPolyPatch::~mappedExtrudedWallPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::mappedExtrudedWallPolyPatch::samplePoints() const
{
    if (!samplePointsPtr_.valid())
    {
        const bool isExtrudedRegion = bottomPatch_ != word::null;

        if (isExtrudedRegion)
        {
            // If this is the extruded region we need to work out where the
            // corresponding sampling points are on the bottom patch. We do
            // this by waving the bottom patch points across the layers.

            const polyMesh& mesh = boundaryMesh().mesh();
            const polyPatch& pp = *this;
            const polyPatch& bottomPp = boundaryMesh()[bottomPatch_];

            // Get the sample points from the bottom patch
            const pointField bottomSamplePoints
            (
                refCast<const mappedPatchBase>(bottomPp).samplePoints()
            );

            // Initialise faces on the bottom patch to wave from
            labelList initialFaces(bottomPp.size());
            List<LayerInfoData<point>> initialFaceInfo(bottomPp.size());
            forAll(bottomPp, bottomPpFacei)
            {
                initialFaces[bottomPpFacei] =
                    bottomPp.start() + bottomPpFacei;
                initialFaceInfo[bottomPpFacei] =
                    LayerInfoData<point>
                    (
                        0,
                        -1,
                        bottomSamplePoints[bottomPpFacei]
                    );
            }

            // Wave across the mesh layers
            List<LayerInfoData<point>> faceInfo(mesh.nFaces());
            List<LayerInfoData<point>> cellInfo(mesh.nCells());
            FaceCellWave<LayerInfoData<point>>
            (
                mesh,
                initialFaces,
                initialFaceInfo,
                faceInfo,
                cellInfo,
                mesh.globalData().nTotalCells() + 1
            );

            // Unpack into this patch's sample points
            samplePointsPtr_.set(new pointField(pp.size()));
            forAll(pp, ppFacei)
            {
                const LayerInfoData<point>& info =
                    faceInfo[pp.start() + ppFacei];

                static nil td;

                if (!info.valid(td))
                {
                    FatalErrorInFunction
                        << "Mesh \"" << mesh.name()
                        << "\" is not layered between the extruded wall patch "
                        << "\"" << pp.name() << "\" and the bottom patch \""
                        << bottomPp.name() << "\"" << exit(FatalError);
                }

                samplePointsPtr_()[ppFacei] = info.data();
            }
        }
        else
        {
            // If this is not the extruded region then we trigger construction
            // of mapping on the extruded region and then reverse map the
            // extruded region's sampling locations so they are available here

            const mappedExtrudedWallPolyPatch& samplePp =
                refCast<const mappedExtrudedWallPolyPatch>(samplePolyPatch());

            samplePointsPtr_.set
            (
                samplePp.reverseDistribute
                (
                    samplePp.mappedPatchBase::samplePoints()
                ).ptr()
            );
        }
    }

    return samplePointsPtr_();
}


void Foam::mappedExtrudedWallPolyPatch::write(Ostream& os) const
{
    mappedWallPolyPatch::write(os);
    writeEntryIfDifferent(os, "bottomPatch", word::null, bottomPatch_);
}


// ************************************************************************* //
