/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "sampledPatch.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace sampledSurfaces
{
    defineTypeNameAndDebug(patch, 0);
    addToRunTimeSelectionTable(sampledSurface, patch, word);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledSurfaces::patch::patch
(
    const word& name,
    const polyMesh& mesh,
    const wordReList& patchNames,
    const bool triangulate
)
:
    sampledSurface(name, mesh),
    patchNames_(patchNames),
    triangulate_(triangulate),
    needsUpdate_(true)
{}


Foam::sampledSurfaces::patch::patch
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledSurface(name, mesh, dict),
    patchNames_(dict.lookup("patches")),
    triangulate_(dict.lookupOrDefault("triangulate", false)),
    needsUpdate_(true)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledSurfaces::patch::~patch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelList& Foam::sampledSurfaces::patch::patchIndices() const
{
    if (patchIndices_.empty())
    {
        patchIndices_ = mesh().boundaryMesh().patchSet
        (
            patchNames_,
            false
        ).sortedToc();
    }
    return patchIndices_;
}


bool Foam::sampledSurfaces::patch::needsUpdate() const
{
    return needsUpdate_;
}


bool Foam::sampledSurfaces::patch::expire()
{
    // already marked as expired
    if (needsUpdate_)
    {
        return false;
    }

    sampledSurface::clearGeom();
    MeshStorage::clear();
    patchIndices_.clear();
    patchIndex_.clear();
    patchFaceLabels_.clear();
    patchStart_.clear();

    needsUpdate_ = true;
    return true;
}


bool Foam::sampledSurfaces::patch::update()
{
    if (!needsUpdate_)
    {
        return false;
    }

    label sz = 0;
    forAll(patchIndices(), i)
    {
        label patchi = patchIndices()[i];
        const polyPatch& pp = mesh().boundaryMesh()[patchi];

        if (isA<emptyPolyPatch>(pp))
        {
            FatalErrorInFunction
                << "Cannot sample an empty patch. Patch " << pp.name()
                << exit(FatalError);
        }

        sz += pp.size();
    }

    // For every face (or triangle) the originating patch and local face in the
    // patch.
    patchIndex_.setSize(sz);
    patchFaceLabels_.setSize(sz);
    patchStart_.setSize(patchIndices().size());
    labelList meshFaceLabels(sz);

    sz = 0;

    forAll(patchIndices(), i)
    {
        label patchi = patchIndices()[i];

        patchStart_[i] = sz;

        const polyPatch& pp = mesh().boundaryMesh()[patchi];

        forAll(pp, j)
        {
            patchIndex_[sz] = i;
            patchFaceLabels_[sz] = j;
            meshFaceLabels[sz] = pp.start()+j;
            sz++;
        }
    }

    indirectPrimitivePatch allPatches
    (
        IndirectList<face>(mesh().faces(), meshFaceLabels),
        mesh().points()
    );

    this->storedPoints() = allPatches.localPoints();
    this->storedFaces()  = allPatches.localFaces();


    // triangulate uses remapFaces()
    // - this is somewhat less efficient since it recopies the faces
    // that we just created, but we probably don't want to do this
    // too often anyhow.
    if (triangulate_)
    {
        MeshStorage::triangulate();
    }

    if (debug)
    {
        print(Pout);
        Pout<< endl;
    }

    needsUpdate_ = false;
    return true;
}


// remap action on triangulation
void Foam::sampledSurfaces::patch::remapFaces(const labelUList& faceMap)
{
    // recalculate the cells cut
    if (notNull(faceMap) && faceMap.size())
    {
        MeshStorage::remapFaces(faceMap);
        patchFaceLabels_ = labelList
        (
            UIndirectList<label>(patchFaceLabels_, faceMap)
        );
        patchIndex_ = labelList
        (
            UIndirectList<label>(patchIndex_, faceMap)
        );

        // Redo patchStart.
        if (patchIndex_.size() > 0)
        {
            patchStart_[patchIndex_[0]] = 0;
            for (label i = 1; i < patchIndex_.size(); i++)
            {
                if (patchIndex_[i] != patchIndex_[i-1])
                {
                    patchStart_[patchIndex_[i]] = i;
                }
            }
        }
    }
}


Foam::tmp<Foam::scalarField> Foam::sampledSurfaces::patch::sample
(
    const volScalarField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::vectorField> Foam::sampledSurfaces::patch::sample
(
    const volVectorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::sphericalTensorField> Foam::sampledSurfaces::patch::sample
(
    const volSphericalTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::symmTensorField> Foam::sampledSurfaces::patch::sample
(
    const volSymmTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::tensorField> Foam::sampledSurfaces::patch::sample
(
    const volTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::scalarField> Foam::sampledSurfaces::patch::sample
(
    const surfaceScalarField& sField
) const
{
    return sampleField(sField);
}


Foam::tmp<Foam::vectorField> Foam::sampledSurfaces::patch::sample
(
    const surfaceVectorField& sField
) const
{
    return sampleField(sField);
}


Foam::tmp<Foam::sphericalTensorField> Foam::sampledSurfaces::patch::sample
(
    const surfaceSphericalTensorField& sField
) const
{
    return sampleField(sField);
}


Foam::tmp<Foam::symmTensorField> Foam::sampledSurfaces::patch::sample
(
    const surfaceSymmTensorField& sField
) const
{
    return sampleField(sField);
}


Foam::tmp<Foam::tensorField> Foam::sampledSurfaces::patch::sample
(
    const surfaceTensorField& sField
) const
{
    return sampleField(sField);
}


Foam::tmp<Foam::scalarField> Foam::sampledSurfaces::patch::interpolate
(
    const interpolation<scalar>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::vectorField> Foam::sampledSurfaces::patch::interpolate
(
    const interpolation<vector>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::sphericalTensorField> Foam::sampledSurfaces::patch::interpolate
(
    const interpolation<sphericalTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::symmTensorField> Foam::sampledSurfaces::patch::interpolate
(
    const interpolation<symmTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::tensorField> Foam::sampledSurfaces::patch::interpolate
(
    const interpolation<tensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


void Foam::sampledSurfaces::patch::print(Ostream& os) const
{
    os  << "patch: " << name() << " :"
        << "  patches:" << patchNames()
        << "  faces:" << faces().size()
        << "  points:" << points().size();
}


// ************************************************************************* //
