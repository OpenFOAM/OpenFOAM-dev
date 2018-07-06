/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::meshDualiser

Description
    Creates dual of polyMesh. Every point becomes a cell (or multiple cells
    for feature points), a walk around every edge creates faces between them.

    Put all points you want in the final mesh into featurePoints; all edge(mid)s
    you want in the final mesh into featureEdges; all face(centre)s in
    faceFaces.

    Usually to preserve boundaries:
        - all boundary faces are featureFaces
        - all edges and points in between different patches are
          featureEdges/points.

    In same way you can also preserve internal faces (e.g. faceZones)

SourceFiles
    Foam::meshDualiser.C

\*---------------------------------------------------------------------------*/

#ifndef meshDualiser_H
#define meshDualiser_H

#include "DynamicList.H"
#include "PackedBoolList.H"
#include "boolList.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;
class polyTopoChange;

/*---------------------------------------------------------------------------*\
                           Class meshDualiser Declaration
\*---------------------------------------------------------------------------*/

class meshDualiser
{
    // Private data

        const polyMesh& mesh_;

        //- From point on cell to dual cell. Either single entry or
        //  one entry per pointCells
        labelListList pointToDualCells_;

        //- From point to dual point (or -1 if not feature point).
        labelList pointToDualPoint_;

        //- From cell to dual point. All cells become point
        labelList cellToDualPoint_;

        //- From face to dual point (or -1 if not feature face)
        labelList faceToDualPoint_;

        //- From edge to dual point (or -1 if not feature edge)
        labelList edgeToDualPoint_;


    // Private Member Functions

        static void checkPolyTopoChange(const polyTopoChange&);

        static void dumpPolyTopoChange(const polyTopoChange&, const fileName&);

        //- Find dual cell given point and cell
        label findDualCell(const label celli, const label pointi) const;

        //- Helper function to generate dualpoints on all boundary edges
        //  emanating from (boundary & feature) point
        void generateDualBoundaryEdges
        (
            const PackedBoolList&,
            const label pointi,
            polyTopoChange&
        );

        //- Check that owner and neighbour of face have same dual cell
        bool sameDualCell
        (
            const label facei,
            const label pointi
        ) const;

        //- Add internal face
        label addInternalFace
        (
            const label masterPointi,
            const label masterEdgeI,
            const label masterFacei,

            const bool edgeOrder,
            const label dualCell0,
            const label dualCell1,
            const DynamicList<label>& verts,
            polyTopoChange& meshMod
        ) const;

        //- Add boundary face
        label addBoundaryFace
        (
            const label masterPointi,
            const label masterEdgeI,
            const label masterFacei,

            const label dualCelli,
            const label patchi,
            const DynamicList<label>& verts,
            polyTopoChange& meshMod
        ) const;

        //- Create internal faces walking around edge
        void createFacesAroundEdge
        (
            const bool splitFace,
            const PackedBoolList&,
            const label edgeI,
            const label startFacei,
            polyTopoChange&,
            boolList& doneEFaces
        ) const;

        //- Create single internal face from internal face
        void createFaceFromInternalFace
        (
            const label facei,
            label& fp,
            polyTopoChange&
        ) const;

        //- Creates boundary faces walking around point on patchi.
        void createFacesAroundBoundaryPoint
        (
            const label patchi,
            const label patchPointi,
            const label startFacei,
            polyTopoChange&,
            boolList& donePFaces            // pFaces visited
        ) const;

        //- Disallow default bitwise copy construct
        meshDualiser(const meshDualiser&);

        //- Disallow default bitwise assignment
        void operator=(const meshDualiser&);


public:

    //- Runtime type information
    ClassName("meshDualiser");


    // Constructors

        //- Construct from mesh
        meshDualiser(const polyMesh&);


    // Member Functions

        // Access

            //- From point on cell to dual cell. Either single entry or
            //  one entry per pointCells.
            const labelListList& pointToDualCells() const
            {
                return pointToDualCells_;
            }

            //- From point to dual point (or -1 if not feature point).
            const labelList& pointToDualPoint() const
            {
                return pointToDualPoint_;
            }

            //- From cell to dual point (at cell centre). All cells become
            //  points.
            const labelList& cellToDualPoint() const
            {
                return cellToDualPoint_;
            }

            //- From face to dual point (at face centre; or -1 if not
            //  feature face).
            const labelList& faceToDualPoint() const
            {
                return faceToDualPoint_;
            }

            //- From edge to dual point (at edge mid; or -1 if not feature
            //  edge).
            const labelList& edgeToDualPoint() const
            {
                return edgeToDualPoint_;
            }


        // Edit


            //- Insert all changes into meshMod to convert the polyMesh into
            //  its dual.
            //  featureFaces  : faces where we want a point at the face centre
            //  featureEdges  : edges           ,,                 edge mid
            //  featurePoints : points          ,,    point. Two variants:
            //      singleCellFeaturePoints : point becomes one dualcell.
            //          Use this for e.g. convex boundary points.
            //      multiCellFeaturePoints : one dualcell per original cell
            //          around point. Use this for e.g. concave boundary points
            //          since it prevents big concave boundary cells.
            void setRefinement
            (
                const bool splitFace,
                const labelList& featureFaces,
                const labelList& featureEdges,
                const labelList& singleCellFeaturePoints,
                const labelList& multiCellFeaturePoints,
                polyTopoChange& meshMod
            );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
