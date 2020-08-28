/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

Class
    Foam::polyDualMesh

Description
    Creates dual of polyMesh.

    Every polyMesh point becomes a cell on the dual mesh

    Every polyMesh cell and patchFace becomes a point on the dual mesh.


SourceFiles
    polyDualMesh.C

\*---------------------------------------------------------------------------*/

#ifndef polyDualMesh_H
#define polyDualMesh_H

#include "polyMesh.H"
#include "labelIOList.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class polyDualMesh Declaration
\*---------------------------------------------------------------------------*/

class polyDualMesh
:
    public polyMesh
{
    // Private Data

        //- From polyMesh cell to my point
        labelIOList cellPoint_;

        //- From polyMesh boundary face (face-mesh.nInternalFaces()) to my point
        labelIOList boundaryFacePoint_;


    // Private Member Functions

        static labelList getFaceOrder
        (
            const labelList& faceOwner,
            const labelList& faceNeighbour,
            const cellList& cells,
            label& nInternalFaces
        );

        static void getPointEdges
        (
            const primitivePatch& patch,
            const label facei,
            const label pointi,
            label& e0,
            label& e1
        );

        static labelList collectPatchSideFace
        (
            const polyPatch& patch,
            const label patchToDualOffset,
            const labelList& edgeToDualPoint,
            const labelList& pointToDualPoint,
            const label pointi,

            label& edgeI
        );

        static void collectPatchInternalFace
        (
            const polyPatch& patch,
            const label patchToDualOffset,
            const labelList& edgeToDualPoint,

            const label pointi,
            const label startEdgeI,

            labelList& dualFace2,
            labelList& featEdgeIndices2
        );

        static void splitFace
        (
            const polyPatch& patch,
            const labelList& pointToDualPoint,

            const label pointi,
            const labelList& dualFace,
            const labelList& featEdgeIndices,

            DynamicList<face>& dualFaces,
            DynamicList<label>& dualOwner,
            DynamicList<label>& dualNeighbour,
            DynamicList<label>& dualRegion
        );

        static void dualPatch
        (
            const polyPatch& patch,
            const label patchToDualOffset,
            const labelList& edgeToDualPoint,
            const labelList& pointToDualPoint,

            const pointField& dualPoints,

            DynamicList<face>& dualFaces,
            DynamicList<label>& dualOwner,
            DynamicList<label>& dualNeighbour,
            DynamicList<label>& dualRegion
        );

        void calcDual
        (
            const polyMesh& mesh,
            const labelList& featureEdges,
            const labelList& featurePoints
        );


public:

    //- Runtime type information
    ClassName("polyDualMesh");


    // Constructors

        //- Construct from IOobject
        polyDualMesh(const IOobject&);

        //- Construct from polyMesh and list of edges and points to represent.
        //  Feature edge and point labels are in local addressing of a patch
        //  over all boundary faces.
        polyDualMesh
        (
            const polyMesh&,
            const labelList& featureEdges,
            const labelList& featurePoints
        );

        //- Construct from polyMesh and feature edge angle. Uses calcFeatures
        //  below to determine feature edges and points.
        polyDualMesh
        (
            const polyMesh&,
            const scalar featureCos
        );

        //- Disallow default bitwise copy construction
        polyDualMesh(const polyDualMesh&) = delete;

        //- Helper function to create feature edges and points based on
        //  feature angle and patches.
        static void calcFeatures
        (
            const polyMesh&,
            const scalar featureCos,
            labelList& featureEdges,
            labelList& featurePoints
        );


    //- Destructor
    ~polyDualMesh();


    // Member Functions

        // Access

            //- From polyMesh cell to dual point
            const labelIOList& cellPoint() const
            {
                return cellPoint_;
            }

            //- From polyMesh patch face to dual point
            const labelIOList& boundaryFacePoint() const
            {
                return boundaryFacePoint_;
            }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const polyDualMesh&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
