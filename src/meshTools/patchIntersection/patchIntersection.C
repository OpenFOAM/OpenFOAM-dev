/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "patchIntersection.H"
#include "primitivePatch.H"
#include "vtkWritePolyData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const bool Foam::patchIntersection::orientToSource_ = true;

namespace Foam
{
    defineTypeNameAndDebug(patchIntersection, 0);
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::patchIntersection::report(const word& writeSuffix)
{
    {
        const primitivePatch patch
        (
            SubList<face>(faces_, faces_.size()),
            points_
        );

        scalar area = 0, srcArea = 0, tgtArea = 0;
        forAll(faces_, facei)
        {
            const scalar a = faces_[facei].mag(points_);
            area += a;
            srcArea += faceSrcFaces_[facei] != -1 ? a : 0;
            tgtArea += faceTgtFaces_[facei] != -1 ? a : 0;
        }
        Info<< indent << "Source/target coverage = " << srcArea/area
            << "/" << tgtArea/area << endl;

        DynamicList<label> nEdgesNFaces, nFacesNEdges;
        forAll(faces_, facei)
        {
            const label n = faces_[facei].size();
            nEdgesNFaces.resize(max(nEdgesNFaces.size(), n + 1), 0);
            ++ nEdgesNFaces[n];
        }
        forAll(patch.edgeFaces(), edgei)
        {
            const label n = patch.edgeFaces()[edgei].size();
            nFacesNEdges.resize(max(nFacesNEdges.size(), n + 1), 0);
            ++ nFacesNEdges[n];
        }
        Info<< indent << "Faces by number of edges = (";
        forAll(nEdgesNFaces, n)
        {
            Info<< (n ? " " : "") << nEdgesNFaces[n];
        }
        Info<< ")" << endl << indent << "Edges by number of faces = (";
        forAll(nFacesNEdges, n)
        {
            Info<< (n ? " " : "") << nFacesNEdges[n];
        }
        Info<< ")" << endl;
    }

    if (debug)
    {
        Info<< indent << "Writing intersected patch" << incrIndent << endl;

        const fileName patchFileName =
            type() + "_patch" + (writeSuffix.empty() ? "" : "_")
          + writeSuffix + ".vtk";
        Info<< indent << "Writing patch to " << patchFileName << endl;
        vtkWritePolyData::write
        (
            patchFileName,
            "intersectedPatch",
            false,
            points_,
            labelList(),
            labelListList(),
            faces_,
            "srcFace",
            false,
            labelField(faceSrcFaces_),
            "tgtFace",
            false,
            labelField(faceTgtFaces_)
        );

        const fileName patchEdgesFileName =
            type() + "_patchEdges" + (writeSuffix.empty() ? "" : "_")
          + writeSuffix + ".vtk";
        Info<< indent << "Writing patch edges to " << patchEdgesFileName
            << endl;
        const primitivePatch patch
        (
            SubList<face>(faces_, faces_.size()),
            points_
        );
        labelField edgeNFaces(patch.nEdges());
        forAll(patch.edgeFaces(), edgei)
        {
            edgeNFaces[edgei] = patch.edgeFaces()[edgei].size();
        }
        vtkWritePolyData::write
        (
            patchEdgesFileName,
            "intersectedPatchEdges",
            false,
            patch.localPoints(),
            labelList(),
            patch.edges(),
            faceList(),
            "nFaces",
            false,
            edgeNFaces
        );

        Info<< decrIndent;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchIntersection::patchIntersection
(
    const label srcNPoints,
    const label tgtNPoints,
    const label srcNEdges,
    const label tgtNEdges,
    const label srcNFaces,
    const label tgtNFaces
)
:
    points_(),

    srcPointPoints_(srcNPoints),
    tgtPointPoints_(tgtNPoints),
    pointSrcPoints_(),
    pointTgtPoints_(),

    srcEdgePoints_(srcNEdges),
    tgtEdgePoints_(tgtNEdges),
    pointSrcEdges_(),
    pointTgtEdges_(),

    pointSrcFaces_(),
    pointTgtFaces_(),

    faces_(),

    srcFaceFaces_(srcNFaces),
    tgtFaceFaces_(tgtNFaces),
    faceSrcFaces_(),
    faceTgtFaces_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::patchIntersection::~patchIntersection()
{}


// ************************************************************************* //
