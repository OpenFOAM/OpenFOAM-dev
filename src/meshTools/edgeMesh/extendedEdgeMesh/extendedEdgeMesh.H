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
    Foam::extendedEdgeMesh

Description

    Description of feature edges and points.

    Feature points are a sorted subset at the start of the overall points list:
        0 .. concaveStart_-1                : convex points (w.r.t normals)
        concaveStart_ .. mixedStart_-1      : concave points
        mixedStart_ .. nonFeatureStart_-1   : mixed internal/external points
        nonFeatureStart_ .. size-1          : non-feature points

    Feature edges are the edgeList of the edgeMesh and are sorted:
        0 .. internalStart_-1           : external edges (convex w.r.t normals)
        internalStart_ .. flatStart_-1  : internal edges (concave)
        flatStart_ .. openStart_-1      : flat edges (neither concave or convex)
                                          can arise from region interfaces on
                                          flat surfaces
        openStart_ .. multipleStart_-1  : open edges (e.g. from baffle surfaces)
        multipleStart_ .. size-1        : multiply connected edges

    The edge direction and feature edge and feature point adjacent normals
    are stored.

SourceFiles
    extendedEdgeMeshI.H
    extendedEdgeMesh.C
    extendedEdgeMeshNew.C

\*---------------------------------------------------------------------------*/

#ifndef extendedEdgeMesh_H
#define extendedEdgeMesh_H

#include "edgeMesh.H"
#include "indexedOctree.H"
#include "treeDataEdge.H"
#include "treeDataPoint.H"
#include "PrimitivePatch.H"
#include "pointIndexHitList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class surfaceFeatures;
class objectRegistry;

// Forward declaration of friend functions and operators

class extendedEdgeMesh;

Istream& operator>>(Istream&, extendedEdgeMesh&);
Ostream& operator<<(Ostream&, const extendedEdgeMesh&);


/*---------------------------------------------------------------------------*\
                      Class extendedEdgeMesh Declaration
\*---------------------------------------------------------------------------*/

class extendedEdgeMesh
:
    public edgeMesh
{

public:

    //- Runtime type information
    TypeName("extendedEdgeMesh");

    enum pointStatus
    {
        CONVEX,     // Fully convex point (w.r.t normals)
        CONCAVE,    // Fully concave point
        MIXED,      // A point surrounded by both convex and concave edges
        NONFEATURE  // Not a feature point
    };

    static const Foam::NamedEnum<pointStatus, 4> pointStatusNames_;

    enum edgeStatus
    {
        EXTERNAL,   // "Convex" edge
        INTERNAL,   // "Concave" edge
        FLAT,       // Neither concave or convex, on a flat surface
        OPEN,       // i.e. only connected to one face
        MULTIPLE,   // Multiply connected (connected to more than two faces)
        NONE        // Not a classified feature edge (consistency with
                    // surfaceFeatures)
    };

    static const Foam::NamedEnum<edgeStatus, 6> edgeStatusNames_;

    //- Normals point to the outside
    enum sideVolumeType
    {
        INSIDE  = 0,  // mesh inside
        OUTSIDE = 1,  // mesh outside
        BOTH    = 2,  // e.g. a baffle
        NEITHER = 3   // not sure when this may be used
    };

    static const Foam::NamedEnum<sideVolumeType, 4> sideVolumeTypeNames_;

    //- Angular closeness tolerance for treating normals as the same
    static scalar cosNormalAngleTol_;


protected:

    // Static data

        //- Index of the start of the convex feature points - static as 0
        static label convexStart_;

        //- Index of the start of the external feature edges - static as 0
        static label externalStart_;


    // Protected data

        //- Index of the start of the concave feature points
        label concaveStart_;

        //- Index of the start of the mixed type feature points
        label mixedStart_;

        //- Index of the start of the non-feature points
        label nonFeatureStart_;

        //- Index of the start of the internal feature edges
        label internalStart_;

        //- Index of the start of the flat feature edges
        label flatStart_;

        //- Index of the start of the open feature edges
        label openStart_;

        //- Index of the start of the multiply-connected feature edges
        label multipleStart_;

        //- Normals of the features, to be referred to by index by both feature
        //  points and edges, unsorted
        vectorField normals_;

        //- Type per normal: which side of normal to mesh
        List<sideVolumeType> normalVolumeTypes_;

        //- Flat and open edges require the direction of the edge
        vectorField edgeDirections_;

        //- Starting directions for the edges.
        //  This vector points to the half of the plane defined by the first
        //  edge normal.
        labelListList normalDirections_;

        //- Indices of the normals that are adjacent to the feature edges
        labelListList edgeNormals_;

        //- Indices of the normals that are adjacent to the feature points
        //  (only valid for 0..nonFeatureStart_-1)
        labelListList featurePointNormals_;

        //- Indices of feature edges attached to feature points. The edges are
        //  ordered so that they can be circulated.
        labelListList featurePointEdges_;

        //- Feature edges which are on the boundary between regions
        labelList regionEdges_;

        //- Search tree for all feature points
        mutable autoPtr<indexedOctree<treeDataPoint>> pointTree_;

        //- Search tree for all edges
        mutable autoPtr<indexedOctree<treeDataEdge>> edgeTree_;

        //- Individual search trees for each type of edge
        mutable PtrList<indexedOctree<treeDataEdge>> edgeTreesByType_;


    // Protected Member Functions

        //- Classify the type of feature point.  Requires valid stored member
        //  data for edges and normals.
        pointStatus classifyFeaturePoint(label ptI) const;

        template<class Patch>
        void sortPointsAndEdges
        (
            const Patch&,
            const labelList& featureEdges,
            const labelList& regionFeatureEdges,
            const labelList& feaurePoints
        );

public:

    // Static data

        //- Number of possible point types (i.e. number of slices)
        static label nPointTypes;

        //- Number of possible feature edge types (i.e. number of slices)
        static label nEdgeTypes;

        //- Can we read this file format?
        static bool canRead(const fileName&, const bool verbose=false);

        //- Can we read this file format?
        static bool canReadType(const word& ext, const bool verbose=false);

        //- Can we write this file format type?
        static bool canWriteType(const word& ext, const bool verbose=false);

        static wordHashSet readTypes();
        static wordHashSet writeTypes();


    // Constructors

        //- Construct null
        extendedEdgeMesh();

        //- Copy constructor
        explicit extendedEdgeMesh(const extendedEdgeMesh&);

        //- Move constructor
        explicit extendedEdgeMesh(extendedEdgeMesh&&);

        //- Construct from file name (uses extension to determine type)
        extendedEdgeMesh(const fileName&);

        //- Construct from file name (uses extension to determine type)
        extendedEdgeMesh(const fileName&, const word& ext);

        //- Construct from Istream
        extendedEdgeMesh(Istream&);

        //- Construct by transferring components (points, edges)
        extendedEdgeMesh
        (
            pointField&&,
            edgeList&&
        );

        //- Construct given a surface with selected edges,points
        //  (surfaceFeatures)
        //  Extracts, classifies and reorders the data from surfaceFeatures.
        extendedEdgeMesh
        (
            const surfaceFeatures& sFeat,
            const boolList& surfBaffleRegions
        );

        //- Construct from PrimitivePatch
        extendedEdgeMesh
        (
            const PrimitivePatch<faceList, pointField>& surf,
            const labelList& featureEdges,
            const labelList& regionFeatureEdges,
            const labelList& featurePoints
        );

        //- Construct from all components
        extendedEdgeMesh
        (
            const pointField& pts,
            const edgeList& eds,
            label concaveStart,
            label mixedStart,
            label nonFeatureStart,
            label internalStart,
            label flatStart,
            label openStart,
            label multipleStart,
            const vectorField& normals,
            const List<sideVolumeType>& normalVolumeTypes,
            const vectorField& edgeDirections,
            const labelListList& normalDirections,
            const labelListList& edgeNormals,
            const labelListList& featurePointNormals,
            const labelListList& featurePointEdges,
            const labelList& regionEdges
        );


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            extendedEdgeMesh,
            fileExtension,
            (
                const fileName& name
            ),
            (name)
        );


    // Selectors

        //- Select constructed from filename (explicit extension)
        static autoPtr<extendedEdgeMesh> New
        (
            const fileName&,
            const word& ext
        );

        //- Select constructed from filename (implicit extension)
        static autoPtr<extendedEdgeMesh> New(const fileName&);


    //- Destructor
    ~extendedEdgeMesh();


    // Member Functions

        // Find

            //- Find nearest surface edge for the sample point.
            void nearestFeaturePoint
            (
                const point& sample,
                scalar searchDistSqr,
                pointIndexHit& info
            ) const;

            //- Find nearest surface edge for the sample point.
            void nearestFeatureEdge
            (
                const point& sample,
                scalar searchDistSqr,
                pointIndexHit& info
            ) const;

            //- Find nearest surface edge for each sample point.
            void nearestFeatureEdge
            (
                const pointField& samples,
                const scalarField& searchDistSqr,
                pointIndexHitList& info
            ) const;

            //- Find the nearest point on each type of feature edge
            void nearestFeatureEdgeByType
            (
                const point& sample,
                const scalarField& searchDistSqr,
                pointIndexHitList& info
            ) const;

            //- Find all the feature points within searchDistSqr of sample
            void allNearestFeaturePoints
            (
                const point& sample,
                scalar searchRadiusSqr,
                pointIndexHitList& info
            ) const;

            //- Find all the feature edges within searchDistSqr of sample
            void allNearestFeatureEdges
            (
                const point& sample,
                const scalar searchRadiusSqr,
                pointIndexHitList& info
            ) const;

            //- Return the minimum distance between disconnected edges
            scalar minDisconnectedDist(const pointIndexHitList& hitList) const;


        // Access

            //- Return the index of the start of the convex feature points
            inline label convexStart() const;

            //- Return the index of the start of the concave feature points
            inline label concaveStart() const;

            //- Return the index of the start of the mixed type feature points
            inline label mixedStart() const;

            //- Return the index of the start of the non-feature points
            inline label nonFeatureStart() const;

            //- Return the index of the start of the external feature edges
            inline label externalStart() const;

            //- Return the index of the start of the internal feature edges
            inline label internalStart() const;

            //- Return the index of the start of the flat feature edges
            inline label flatStart() const;

            //- Return the index of the start of the open feature edges
            inline label openStart() const;

            //- Return the index of the start of the multiply-connected feature
            //  edges
            inline label multipleStart() const;

            //- Return whether or not the point index is a feature point
            inline bool featurePoint(label ptI) const;

            //- Return the normals of the surfaces adjacent to the feature edges
            //  and points
            inline const vectorField& normals() const;

            //- Return
            inline const List<sideVolumeType>& normalVolumeTypes() const;

            //- Return the edgeDirection vectors
            inline const vectorField& edgeDirections() const;

            //-
            inline const labelListList& normalDirections() const;

            //- Return the direction of edgeI, pointing away from ptI
            inline vector edgeDirection(label edgeI, label ptI) const;

            //- Return the indices of the normals that are adjacent to the
            //  feature edges
            inline const labelListList& edgeNormals() const;

            //- Return the normal vectors for a given set of normal indices
            inline vectorField edgeNormals(const labelList& edgeNormIs) const;

            //- Return the normal vectors for a given edge
            inline vectorField edgeNormals(label edgeI) const;

            //- Return the indices of the normals that are adjacent to the
            //  feature points
            inline const labelListList& featurePointNormals() const;

            //- Return the normal vectors for a given feature point
            inline vectorField featurePointNormals(label ptI) const;

            //- Return the edge labels for a given feature point. Edges are
            //  ordered by the faces that they share. The edge labels
            //  correspond to the entry in edges().
            inline const labelListList& featurePointEdges() const;

            //- Return the feature edges which are on the boundary between
            //  regions
            inline const labelList& regionEdges() const;

            //- Return the pointStatus of a specified point
            inline pointStatus getPointStatus(label ptI) const;

            //- Return the edgeStatus of a specified edge
            inline edgeStatus getEdgeStatus(label edgeI) const;

            //- Return the baffle faces of a specified edge
            inline PackedList<2> edgeBaffles(label edgeI) const;

            //- Demand driven construction of octree for feature points
            const indexedOctree<treeDataPoint>& pointTree() const;

            //- Demand driven construction of octree for boundary edges
            const indexedOctree<treeDataEdge>& edgeTree() const;

            //- Demand driven construction of octree for boundary edges by type
            const PtrList<indexedOctree<treeDataEdge>>&
            edgeTreesByType() const;


        // Edit

            //- Transfer the contents of the argument and annul the argument
            void transfer(extendedEdgeMesh&);

            //- Clear all storage
            virtual void clear();

            //- Add extendedEdgeMesh. No filtering of duplicates.
            void add(const extendedEdgeMesh&);

            //- Flip normals. All concave become convex, all internal external
            //  etc.
            void flipNormals();


        // Read

            //- Read from file. Chooses reader based on explicit extension
            bool read(const fileName&, const word& ext);

            //- Read from file. Chooses reader based on detected extension
            virtual bool read(const fileName&);


        // Write

            //- Write all components of the extendedEdgeMesh as obj files
            void writeObj
            (
                const fileName& prefix,
                const bool verbose = true
            ) const;

            //- Dump some information
            virtual void writeStats(Ostream& os) const;

            friend Istream& operator>>(Istream& is, sideVolumeType& vt);
            friend Ostream& operator<<(Ostream& os, const sideVolumeType& vt);


        //- Classify the type of feature edge.  Requires face centre 0 to face
        //  centre 1 vector to distinguish internal from external
        static edgeStatus classifyEdge
        (
            const List<vector>& norms,
            const labelList& edNorms,
            const vector& fC0tofC1
        );


        // IOstream Operators

            friend Ostream& operator<<(Ostream&, const extendedEdgeMesh&);
            friend Istream& operator>>(Istream&, extendedEdgeMesh&);
};


Istream& operator>>(Istream& is, extendedEdgeMesh::sideVolumeType& vt);
Ostream& operator<<(Ostream& os, const extendedEdgeMesh::sideVolumeType& vt);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "extendedEdgeMeshI.H"

#ifdef NoRepository
    #include "extendedEdgeMeshTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
