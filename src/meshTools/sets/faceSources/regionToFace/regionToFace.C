/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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
#include "mappedPatchBase.H"
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

addToRunTimeSelectionTable(topoSetSource, regionToFace, istream);

}


Foam::topoSetSource::addToUsageTable Foam::regionToFace::usage_
(
    regionToFace::typeName,
    "\n    Usage: regionToFace <faceSet> (x y z)\n\n"
    "    Select all faces in the connected region of the faceSet"
    " starting from the point.\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::regionToFace::markZone
(
    const indirectPrimitivePatch& patch,
    const label procI,
    const label faceI,
    const label zoneI,
    labelList& faceZone
) const
{
    // Data on all edges and faces
    List<patchEdgeFaceRegion> allEdgeInfo(patch.nEdges());
    List<patchEdgeFaceRegion> allFaceInfo(patch.size());

    DynamicList<label> changedEdges;
    DynamicList<patchEdgeFaceRegion> changedInfo;

    if (Pstream::myProcNo() == procI)
    {
        const labelList& fEdges = patch.faceEdges()[faceI];
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

    forAll(allFaceInfo, faceI)
    {
        if (allFaceInfo[faceI].region() == zoneI)
        {
            faceZone[faceI] = zoneI;
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

    mappedPatchBase::nearInfo ni
    (
        pointIndexHit(false, vector::zero, -1),
        Tuple2<scalar, label>
        (
            sqr(GREAT),
            Pstream::myProcNo()
        )
    );

    forAll(patch, i)
    {
        const point& fc = patch.faceCentres()[i];
        scalar d2 = magSqr(fc-nearPoint_);

        if (!ni.first().hit() || d2 < ni.second().first())
        {
            ni.second().first() = d2;
            ni.first().setHit();
            ni.first().setPoint(fc);
            ni.first().setIndex(i);
        }
    }

    // Globally reduce
    combineReduce(ni, mappedPatchBase::nearestEqOp());

    Info<< "    Found nearest face at " << ni.first().rawPoint()
        << " on processor " << ni.second().second()
        << " face " << ni.first().index()
        << " distance " << Foam::sqrt(ni.second().first()) << endl;

    labelList faceRegion(patch.size(), -1);
    markZone
    (
        patch,
        ni.second().second(),   // procI
        ni.first().index(),     // start face
        0,                      // currentZone
        faceRegion
    );

    forAll(faceRegion, faceI)
    {
        if (faceRegion[faceI] == 0)
        {
            addOrDelete(set, patch.addressing()[faceI], add);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
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


// Construct from dictionary
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


// Construct from Istream
Foam::regionToFace::regionToFace
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetSource(mesh),
    setName_(checkIs(is)),
    nearPoint_(checkIs(is))
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
