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

#include "mappedExtrudedPatchBase.H"
#include "LayerInfoData.H"
#include "FaceCellWave.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mappedExtrudedPatchBase, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::vectorField>
Foam::mappedExtrudedPatchBase::patchFaceAreas() const
{
    if (!bottomFaceAreasPtr_.valid())
    {
        const bool isExtrudedRegion = oppositePatch_ != word::null;

        if (isExtrudedRegion)
        {
            // If this is the extruded region we need to work out what the
            // corresponding areas and centres are on the bottom patch. We do
            // this by waving these values across the layers.

            const polyMesh& mesh = patch_.boundaryMesh().mesh();
            const polyPatch& pp = patch_;
            const polyPatch& bottomPp = patch_.boundaryMesh()[oppositePatch_];

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
                        << "\" is not layered between the extruded patch "
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

            const mappedExtrudedPatchBase& nbrPp =
                refCast<const mappedExtrudedPatchBase>(nbrPolyPatch());

            bottomFaceAreasPtr_.set
            (
                nbrPp.reverseDistribute
                (
                    nbrPp.patch_.primitivePatch::faceAreas()
                ).ptr()
            );
            bottomFaceCentresPtr_.set
            (
                nbrPp.reverseDistribute
                (
                    nbrPp.patch_.primitivePatch::faceCentres()
                ).ptr()
            );
        }
    }

    return bottomFaceAreasPtr_();
}


Foam::tmp<Foam::pointField>
Foam::mappedExtrudedPatchBase::patchFaceCentres() const
{
    if (!bottomFaceCentresPtr_.valid())
    {
        patchFaceAreas();
    }

    return bottomFaceCentresPtr_();
}


Foam::tmp<Foam::pointField>
Foam::mappedExtrudedPatchBase::patchLocalPoints() const
{
    NotImplemented;
    return tmp<pointField>(nullptr);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mappedExtrudedPatchBase::mappedExtrudedPatchBase(const polyPatch& pp)
:
    mappedPatchBase(pp),
    oppositePatch_(word::null)
{}


Foam::mappedExtrudedPatchBase::mappedExtrudedPatchBase
(
    const polyPatch& pp,
    const word& nbrRegionName,
    const word& nbrPatchName,
    const word& oppositePatch,
    const cyclicTransform& transform
)
:
    mappedPatchBase(pp, nbrRegionName, nbrPatchName, transform),
    oppositePatch_(oppositePatch)
{}


Foam::mappedExtrudedPatchBase::mappedExtrudedPatchBase
(
    const polyPatch& pp,
    const dictionary& dict,
    const bool transformIsNone
)
:
    mappedPatchBase(pp, dict, transformIsNone),
    oppositePatch_(dict.lookupOrDefault<word>("oppositePatch", word::null))
{}


Foam::mappedExtrudedPatchBase::mappedExtrudedPatchBase
(
    const polyPatch& pp,
    const mappedExtrudedPatchBase& mepb
)
:
    mappedPatchBase(pp),
    oppositePatch_(mepb.oppositePatch_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mappedExtrudedPatchBase::~mappedExtrudedPatchBase()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mappedExtrudedPatchBase::clearOut()
{
    mappedPatchBase::clearOut();
    if (reMapAfterMove_)
    {
        bottomFaceCentresPtr_.clear();
    }
}


void Foam::mappedExtrudedPatchBase::write(Ostream& os) const
{
    mappedPatchBase::write(os);
    writeEntryIfDifferent(os, "oppositePatch", word::null, oppositePatch_);
}


// ************************************************************************* //
