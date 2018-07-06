/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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
    Foam::featurePointConformer

Description
    The Delaunay vertices required to conform to a feature point can be
    determined upon initialisation because the feature points are fixed and
    do not change throughout the meshing process.

SourceFiles
    featurePointConformerI.H
    featurePointConformer.C
    featurePointConformerSpecialisations.C

\*---------------------------------------------------------------------------*/

#ifndef featurePointConformer_H
#define featurePointConformer_H

#include "CGALTriangulation3Ddefs.H"
#include "vector.H"
#include "DynamicList.H"
#include "List.H"
#include "extendedFeatureEdgeMesh.H"
#include "pointPairs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class conformalVoronoiMesh;
class cvControls;
class conformationSurfaces;
class pointFeatureEdgesTypes;
class backgroundMeshDecomposition;

/*---------------------------------------------------------------------------*\
                   Class featurePointConformer Declaration
\*---------------------------------------------------------------------------*/

class featurePointConformer
{
    // Static data

        //- Tolerance within which two lines are said to be parallel.
        static const scalar tolParallel;


    // Private data

        //- Reference to the mesher.
        const conformalVoronoiMesh& foamyHexMesh_;

        //- Reference to the mesher controls.
        const cvControls& foamyHexMeshControls_;

        //- Reference to the conformation surfaces.
        const conformationSurfaces& geometryToConformTo_;

        //- Store the feature constraining points, to be reinserted after a
        //  triangulation clear.
        List<Vb> featurePointVertices_;

        //-
        mutable pointPairs<Delaunay> ftPtPairs_;


    // Private Member Functions

        //- Calculate the shared face normal between two edges geometrically.
        vector sharedFaceNormal
        (
            const extendedFeatureEdgeMesh& feMesh,
            const label edgeI,
            const label nextEdgeI
        ) const;

        label getSign(const extendedFeatureEdgeMesh::edgeStatus eStatus) const;

        bool createSpecialisedFeaturePoint
        (
            const extendedFeatureEdgeMesh& feMesh,
            const labelList& pEds,
            const pointFeatureEdgesTypes& pFEdgesTypes,
            const List<extendedFeatureEdgeMesh::edgeStatus>& allEdStat,
            const label ptI,
            DynamicList<Vb>& pts
        ) const;


        void addMasterAndSlavePoints
        (
            const DynamicList<point>& masterPoints,
            const DynamicList<indexedVertexEnum::vertexType>& masterPointsTypes,
            const Map<DynamicList<autoPtr<plane>>>& masterPointReflections,
            DynamicList<Vb>& pts,
            const label ptI
        ) const;

        //- Helper function for conforming to feature points
        void createMasterAndSlavePoints
        (
            const extendedFeatureEdgeMesh& feMesh,
            const label ptI,
            DynamicList<Vb>& pts
        ) const;

        void createMixedFeaturePoints(DynamicList<Vb>& pts) const;

        //- Create the points that will conform to the feature
        void createFeaturePoints(DynamicList<Vb>& pts);

        //- Disallow default bitwise copy construct
        featurePointConformer(const featurePointConformer&);

        //- Disallow default bitwise assignment
        void operator=(const featurePointConformer&);


public:

    //- Runtime type information
    ClassName("featurePointConformer");


    // Constructors

        //- Construct from components
        explicit featurePointConformer
        (
            const conformalVoronoiMesh& foamyHexMesh
        );


    //- Destructor
    ~featurePointConformer();


    // Member Functions

        // Access

            //- Return the feature point vertices for insertion into the
            //  triangulation.
            inline const List<Vb>& featurePointVertices() const;

            //- Return the feature point pair table
            inline const pointPairs<Delaunay>& featurePointPairs() const;


        // Edit

            //- Distribute the feature point vertices according to the
            //  supplied background mesh
            void distribute(const backgroundMeshDecomposition& decomposition);

            //- Reindex the feature point pairs using the map.
            void reIndexPointPairs(const Map<label>& oldToNewIndices);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "featurePointConformerI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
