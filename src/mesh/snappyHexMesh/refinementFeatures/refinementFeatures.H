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
    Foam::refinementFeatures

Description
    Encapsulates queries for features.

SourceFiles
    refinementFeatures.C

\*---------------------------------------------------------------------------*/

#ifndef refinementFeatures_H
#define refinementFeatures_H

#include "extendedFeatureEdgeMesh.H"
#include "indexedOctree.H"
#include "treeDataEdge.H"
#include "treeDataPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class refinementFeatures Declaration
\*---------------------------------------------------------------------------*/

class refinementFeatures
:
    public PtrList<extendedFeatureEdgeMesh>
{
    // Private data

        //- Per shell the list of ranges
        List<scalarField> distances_;

        //- Per shell per distance the refinement level
        labelListList levels_;

        //- Edge
        PtrList<indexedOctree<treeDataEdge>> edgeTrees_;

        //- Features points
        PtrList<indexedOctree<treeDataPoint>> pointTrees_;

        //- Region edge trees (demand driven)
        mutable autoPtr<PtrList<indexedOctree<treeDataEdge>>>
            regionEdgeTreesPtr_;


    // Private Member Functions

        //- Read set of feature edge meshes
        void read(const objectRegistry&, const PtrList<dictionary>&);

        //- Build edge tree and feature point tree
        void buildTrees(const label);

        //- Find shell level higher than ptLevel
        void findHigherLevel
        (
            const pointField& pt,
            const label featI,
            labelList& maxLevel
        ) const;


protected:

        const PtrList<indexedOctree<treeDataEdge>>& edgeTrees() const
        {
            return edgeTrees_;
        }

        const PtrList<indexedOctree<treeDataPoint>>& pointTrees() const
        {
            return pointTrees_;
        }


        const PtrList<indexedOctree<treeDataEdge>>& regionEdgeTrees() const;

public:

    // Constructors

        //- Construct from description
        refinementFeatures
        (
            const objectRegistry& io,
            const PtrList<dictionary>& featDicts
        );


    // Member Functions

        // Access

            //- Per featureEdgeMesh the list of level
            const labelListList& levels() const
            {
                return levels_;
            }

            //- Per featureEdgeMesh the list of ranges
            const List<scalarField>& distances() const
            {
                return distances_;
            }


        // Query

            //- Highest distance of all features
            scalar maxDistance() const;

            //- Find nearest point on nearest feature edge. Sets:
            //
            //    - nearFeature: index of feature mesh
            //    - nearInfo   : location on feature edge and edge index
            //        (note: not feature edge index but index into edges()
            //        directly)
            //    - nearNormal : local feature edge normal
            void findNearestEdge
            (
                const pointField& samples,
                const scalarField& nearestDistSqr,
                labelList& nearFeature,
                List<pointIndexHit>& nearInfo,
                vectorField& nearNormal
            ) const;

            //- Find nearest point on nearest region edge. Sets:
            //
            //    - nearFeature: index of feature mesh
            //    - nearInfo   : location on feature edge and edge index
            //        (note: not feature edge index but index into edges()
            //        directly)
            //    - nearNormal : local feature edge normal
            void findNearestRegionEdge
            (
                const pointField& samples,
                const scalarField& nearestDistSqr,
                labelList& nearFeature,
                List<pointIndexHit>& nearInfo,
                vectorField& nearNormal
            ) const;

            //- Find nearest feature point. Sets:
            //
            //    - nearFeature: index of feature mesh
            //    - nearInfo   : location on feature point and point index.
            //        (note: not index into shapes().pointLabels() but index
            //        into points() directly)
            void findNearestPoint
            (
                const pointField& samples,
                const scalarField& nearestDistSqr,
                labelList& nearFeature,
                List<pointIndexHit>& nearInfo
            ) const;

            //- Find shell level higher than ptLevel
            void findHigherLevel
            (
                const pointField& pt,
                const labelList& ptLevel,
                labelList& maxLevel
            ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
