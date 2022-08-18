/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2022 OpenFOAM Foundation
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

#include "regionToFace.H"
#include "polyMesh.H"
#include "faceSet.H"
#include "RemoteData.H"
#include "indirectPrimitivePatch.H"
#include "PatchTools.H"
#include "addToRunTimeSelectionTable.H"
#include "PatchEdgeFaceWave.H"
#include "patchEdgeFaceRegion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(regionToFace, 0);
    addToRunTimeSelectionTable(topoSetSource, regionToFace, word);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::regionToFace::markZone
(
    const indirectPrimitivePatch& patch,
    const label proci,
    const label facei,
    const label zoneI,
    labelList& faceZone
) const
{
    // Data on all edges and faces
    List<patchEdgeFaceRegion> allEdgeInfo(patch.nEdges());
    List<patchEdgeFaceRegion> allFaceInfo(patch.size());

    DynamicList<label> changedEdges;
    DynamicList<patchEdgeFaceRegion> changedInfo;

    if (Pstream::myProcNo() == proci)
    {
        const labelList& fEdges = patch.faceEdges()[facei];
        forAll(fEdges, i)
        {
            changedEdges.append(fEdges[i]);
            changedInfo.append(zoneI);
        }
    }

    // Walk
    PatchEdgeFaceWave
    <
        indirectPrimitivePatch,
        patchEdgeFaceRegion
    > calc
    (
        mesh_,
        patch,
        changedEdges,
        changedInfo,
        allEdgeInfo,
        allFaceInfo,
        returnReduce(patch.nEdges(), sumOp<label>())
    );

    forAll(allFaceInfo, facei)
    {
        if (allFaceInfo[facei].region() == zoneI)
        {
            faceZone[facei] = zoneI;
        }
    }
}


void Foam::regionToFace::combine(topoSet& set, const bool add) const
{
    Info<< "    Loading subset " << setName_ << " to delimit search region."
        << endl;
    faceSet subSet(mesh_, setName_);

    indirectPrimitivePatch patch
    (
        IndirectList<face>(mesh_.faces(), subSet.toc()),
        mesh_.points()
    );

    RemoteData<scalar> ni;

    forAll(patch, i)
    {
        const point& fc = patch.faceCentres()[i];

        const scalar dSqr = magSqr(fc - nearPoint_);

        if (ni.proci == -1 || dSqr < ni.data)
        {
            ni.proci = Pstream::myProcNo();
            ni.elementi = i;
            ni.data = dSqr;
        }
    }

    // Globally reduce
    combineReduce(ni, RemoteData<scalar>::smallestEqOp());

    Info<< "    Found nearest face on processor " << ni.proci
        << " face " << ni.elementi
        << " distance " << Foam::sqrt(ni.data) << endl;

    labelList faceRegion(patch.size(), -1);
    markZone
    (
        patch,
        ni.proci,
        ni.elementi,
        0,
        faceRegion
    );

    forAll(faceRegion, facei)
    {
        if (faceRegion[facei] == 0)
        {
            addOrDelete(set, patch.addressing()[facei], add);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionToFace::regionToFace
(
    const polyMesh& mesh,
    const word& setName,
    const point& nearPoint
)
:
    topoSetSource(mesh),
    setName_(setName),
    nearPoint_(nearPoint)
{}


Foam::regionToFace::regionToFace
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    setName_(dict.lookup("set")),
    nearPoint_(dict.lookup("nearPoint"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionToFace::~regionToFace()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionToFace::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding all faces of connected region of set "
            << setName_
            << " starting from point "
            << nearPoint_ << " ..." << endl;

        combine(set, true);
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing all cells of connected region of set "
            << setName_
            << " starting from point "
            << nearPoint_ << " ..." << endl;

        combine(set, false);
    }
}


// ************************************************************************* //
