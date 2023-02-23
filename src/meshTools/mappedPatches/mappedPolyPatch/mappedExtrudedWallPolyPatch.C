/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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
    bottomFaceCentresPtr_.clear();
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
    bottomFaceCentresPtr_.clear();
}


void Foam::mappedExtrudedWallPolyPatch::initTopoChange(PstreamBuffers& pBufs)
{
    mappedWallPolyPatch::initTopoChange(pBufs);
}


void Foam::mappedExtrudedWallPolyPatch::topoChange(PstreamBuffers& pBufs)
{
    mappedWallPolyPatch::topoChange(pBufs);
    bottomFaceCentresPtr_.clear();
}


Foam::tmp<Foam::vectorField>
Foam::mappedExtrudedWallPolyPatch::patchFaceAreas() const
{
    if (!bottomFaceAreasPtr_.valid())
    {
        const bool isExtrudedRegion = bottomPatch_ != word::null;

        if (isExtrudedRegion)
        {
            // If this is the extruded region we need to work out what the
            // corresponding areas and centres are on the bottom patch. We do
            // this by waving these values across the layers.

            const polyMesh& mesh = boundaryMesh().mesh();
            const polyPatch& pp = *this;
            const polyPatch& bottomPp = boundaryMesh()[bottomPatch_];

            // Initialise faces on the bottom patch to wave from
            labelList initialFaces(bottomPp.size());
            List<LayerInfoData<Pair<vector>>> initialFaceInfo(bottomPp.size());
            forAll(bottomPp, bottomPpFacei)
            {
                initialFaces[bottomPpFacei] =
                    bottomPp.start() + bottomPpFacei;
                initialFaceInfo[bottomPpFacei] =
                    LayerInfoData<Pair<vector>>
                    (
                        0,
                        -1,
                        Pair<vector>
                        (
                            bottomPp.faceAreas()[bottomPpFacei],
                            bottomPp.faceCentres()[bottomPpFacei]
                        )
                    );
            }

            // Wave across the mesh layers
            List<LayerInfoData<Pair<vector>>> faceInfo(mesh.nFaces());
            List<LayerInfoData<Pair<vector>>> cellInfo(mesh.nCells());
            FaceCellWave<LayerInfoData<Pair<vector>>>
            (
                mesh,
                initialFaces,
                initialFaceInfo,
                faceInfo,
                cellInfo,
                mesh.globalData().nTotalCells() + 1
            );

            // Unpack into this patch's bottom face areas and centres
            bottomFaceAreasPtr_.set(new vectorField(pp.size()));
            bottomFaceCentresPtr_.set(new pointField(pp.size()));
            forAll(pp, ppFacei)
            {
                const LayerInfoData<Pair<vector>>& info =
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

                bottomFaceAreasPtr_()[ppFacei] = info.data().first();
                bottomFaceCentresPtr_()[ppFacei] = info.data().second();
            }
        }
        else
        {
            // If this is not the extruded region then we trigger construction
            // of mapping on the extruded region and then reverse map the
            // extruded region's data so it is available here

            const mappedExtrudedWallPolyPatch& nbrPp =
                refCast<const mappedExtrudedWallPolyPatch>(nbrPolyPatch());

            bottomFaceAreasPtr_.set
            (
                nbrPp.reverseDistribute
                (
                    nbrPp.primitivePatch::faceAreas()
                ).ptr()
            );
            bottomFaceCentresPtr_.set
            (
                nbrPp.reverseDistribute
                (
                    nbrPp.primitivePatch::faceCentres()
                ).ptr()
            );
        }
    }

    return bottomFaceAreasPtr_();
}


Foam::tmp<Foam::pointField>
Foam::mappedExtrudedWallPolyPatch::patchFaceCentres() const
{
    if (!bottomFaceCentresPtr_.valid())
    {
        patchFaceAreas();
    }

    return bottomFaceCentresPtr_();
}


Foam::tmp<Foam::pointField>
Foam::mappedExtrudedWallPolyPatch::patchLocalPoints() const
{
    NotImplemented;
    return tmp<pointField>(nullptr);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

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

void Foam::mappedExtrudedWallPolyPatch::write(Ostream& os) const
{
    mappedWallPolyPatch::write(os);
    writeEntryIfDifferent(os, "bottomPatch", word::null, bottomPatch_);
}


// ************************************************************************* //
