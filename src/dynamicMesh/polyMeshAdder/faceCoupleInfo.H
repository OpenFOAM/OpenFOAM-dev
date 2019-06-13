/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::faceCoupleInfo

Description
    Container for information needed to couple to meshes. When constructed
    from two meshes and a geometric tolerance finds the corresponding
    boundary faces.

    The information it keeps is the set of faces&points (cutFaces,
    cutPoints) that should replace a set of faces on the master
    (masterPatch) and a set of faces on the slave (slavePatch)


    Uses same tolerance to match faces and points on matched faces since
    they both originate from the same points and the tolerance usually
    comes from writing these points with limited precision (6 by default)

    -# Perfect match:
       - one-to-one match for faces and points.
       - the cut is always the 'most connected' of the master and slave so
         multiple master or slave points might point to the same cut point.

    \verbatim
    e.g. master:

        +--+
        |  |
        |  |
        +--+
            +--+
            |  |
            |  |
            +--+
    slave:
            +--+
            |  |
            |  |
            +--+
        +--+
        |  |
        |  |
        +--+
    \endverbatim
    adding both together creates a singly connected 2x2 cavity so suddenly
    the duplicate master points and the duplicate slave points all become
    a single cut point.


    -# Subdivision match:
       - Can be constructed from slave being subdivision of master with the
         polyPatch constructor.
       - Does not include above shared-point detection!

    Notes on multiple slave faces per master:

    As long as
    - all master edges are present in slave
    - slave can have extra edges/points/faces BUT all subfaces have to have
      at least one point on a maste face.

    \verbatim
    So master:
    +-------+
    |       |
    |       |
    |       |
    |       |
    |       |
    |       |
    |       |
    +-------+

    slave:
    +---+---+
    |\  |  /|
    | \ | / |
    |  \|/  |
    +---+---+
    |  /|\  |
    | / | \ |
    |/  |  \|
    +---+---+
    is ok.
    \endverbatim

    For this kind of matching the order is:
    - match cutpoint to masterpoint
    - find those cutEdges that align with a master edge. This gives two sets
      of cut edges: those that have a master equivalent ('border edges') and
      those that don't ('internal edges'). The border edges now divide the
      cutFaces into regions with the same masterFace correspondence.
    - find cutFaces that are fully determined by the border edges they use.
    - all cutFaces that are connected through an internal edge have the same
      master face.


    Note: matching refined faces onto master is a bit dodgy and will probably
    only work for unwarped faces. Also it will fail if e.g. face is split
    into 3x3 since then middle face has no point/edge in common with master.
    (problem is in face matching (findSlavesCoveringMaster), probably
     point/edge matching might just work)


SourceFiles
    faceCoupleInfo.C


\*---------------------------------------------------------------------------*/

#ifndef faceCoupleInfo_H
#define faceCoupleInfo_H

#include "pointField.H"
#include "indirectPrimitivePatch.H"
#include "primitiveFacePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

typedef HashTable<labelList, edge, Hash<edge>> edgeLookup;


// Forward declaration of classes
class face;
class primitiveMesh;
class polyPatch;
class polyMesh;

/*---------------------------------------------------------------------------*\
                           Class faceCoupleInfo Declaration
\*---------------------------------------------------------------------------*/

class faceCoupleInfo
{
    // Private Data

        //- Angle matching tolerance.
        static const scalar angleTol_;

        //- Master patch
        autoPtr<indirectPrimitivePatch> masterPatchPtr_;

        //- Slave patch
        autoPtr<indirectPrimitivePatch> slavePatchPtr_;


        //- Description of cut.
        //  - Cut is the matching area between the slave
        //  and the master.
        //  - cut is the finest of master and slave. It can never be
        //  coarser than either one of them. (so face addressing we keep is
        //  cut-to-master and cut-to-slave)
        //  - multiple master or slave points can end up becoming one cut point
        //  (so point addressing we keep is master-to-cut and slave-to-cut)

        //  Cut consists of faces and points (note: could be expressed as some
        //  kind of PrimitivePatch which holds points instead of reference to
        //  them)
        //  Orientation of cutFaces should be same as masterFaces!
        pointField cutPoints_;
        autoPtr<primitiveFacePatch> cutFacesPtr_;

        //- Additional point coupling information. Is between points on
        //  boundary of both meshes.

        // Addressing to/from cut

            //- master
            labelList cutToMasterFaces_;
            labelList masterToCutPoints_;

            //- slave
            labelList cutToSlaveFaces_;
            labelList slaveToCutPoints_;

            //- For edges originating from splitting of edges:
            //  given the two endpoints of the unsplit edge give the list
            //  of in between vertices
            edgeLookup cutEdgeToPoints_;


    // Private Member Functions

        // Debugging

            //- Calculate face centres from (subset of) faces.
            template<template<class> class FaceList>
            static pointField calcFaceCentres
            (
                const FaceList<face>&,
                const pointField&,
                const label start,
                const label size
            );

            //- Calculate face point averages from (subset of) faces.
            template<template<class> class FaceList>
            static pointField calcFacePointAverages
            (
                const FaceList<face>&,
                const pointField&,
                const label start,
                const label size
            );

            //- Write edges
            static void writeOBJ
            (
                const fileName& fName,
                const edgeList& edges,
                const pointField& points,
                const bool compact = true
            );

            //- Write edges
            static void writeOBJ
            (
                const fileName& fName,
                const pointField& points0,
                const pointField& points1
            );

            //- Write connections between corresponding points and faces
            //  as .obj files.
            void writePointsFaces() const;

            //- Write connections between corresponding edges as .obj files.
            void writeEdges(const labelList&, const labelList&) const;


        // Edge handling/matching

            //- Find corresponding edges on patch when having only a map for
            //  the points.
            labelList findMappedEdges
            (
                const edgeList& edges,
                const labelList& pointMap,
                const indirectPrimitivePatch&
            );

            //- Check if edge on slavePatch corresponds to an edge between faces
            //  in two different polyPatches on the mesh.
            bool regionEdge(const polyMesh&, const label slaveEdgeI) const;

            //- Finds edge connected to point most aligned with master edge.
            label mostAlignedCutEdge
            (
                const bool report,
                const polyMesh& slaveMesh,
                const bool patchDivision,
                const labelList& cutToMasterEdges,
                const labelList& cutToSlaveEdges,
                const label pointi,
                const label edgeStart,
                const label edgeEnd
            ) const;

            //- From (many-to-one) map of cut edges to master edges determine
            //  points in between. I.e. just string up the edges. Stores this
            //  all on cutEdgeToPoints_
            void setCutEdgeToPoints(const labelList& cutToMasterEdges);

        // Face matching

            //- Matches two faces.
            //  Determines rotation for f1 to match up with f0,
            //  i.e. the index in f0 of the first point of f1.
            static label matchFaces
            (
                const scalar absTol,
                const pointField& points0,
                const face& f0,
                const pointField& points1,
                const face& f1,
                const bool sameOrientation
            );

            //- Matches points on patch to points on cut.
            static bool matchPointsThroughFaces
            (
                const scalar absTol,
                const pointField& cutPoints,
                const faceList& cutFaces,
                const pointField& patchPoints,
                const faceList& patchFaces,
                const bool sameOrientation,

                labelList& patchToCutPoints,// patch to (uncompacted) cut points
                labelList& cutToCompact,    // compaction list
                labelList& compactToCut     // compaction list
            );

            //- Returns max distance to masterF of any point on cutF.
            static scalar maxDistance
            (
                const face& cutF,
                const pointField& cutPoints,
                const face& masterF,
                const pointField& masterPoints
            );

            //- Finds matching (boundary)face centres.
            //  Since faces identical uses geometric match on face centres.
            static void findPerfectMatchingFaces
            (
                const primitiveMesh& mesh0,
                const primitiveMesh& mesh1,
                const scalar absTol,

                labelList& mesh0Faces,
                labelList& mesh1Faces
            );

            //- Find matching (boundary)faces. Matching if slave is on top of
            //  master face (slaves is subdivision of master)
            static void findSlavesCoveringMaster
            (
                const primitiveMesh& mesh0,
                const primitiveMesh& mesh1,
                const scalar absTol,

                labelList& mesh0Faces,
                labelList& mesh1Faces
            );

            //- Grow cutToMasterFace across 'internal' edges.
            label growCutFaces(const labelList&, Map<labelList>&);

            void checkMatch(const labelList& cutToMasterEdges) const;

            //- Gets a list of cutFaces (that use a master edge) and the
            // candidate master faces.
            // Checks among these master faces if there is only one remaining
            // unmatched one.
            label matchEdgeFaces(const labelList&, Map<labelList>& candidates);

            //- Gets a list of cutFaces (that use a master edge) and the
            //  candidate master faces.
            //  Finds most aligned master face.
            label geometricMatchEdgeFaces(Map<labelList>& candidates);

        //- Used by perfectPointMatch. Determine match from cut points to
        //  slave points (for perfect matching faces)
        void perfectSlavePointMatch(const scalar absTol);

        //- Find point and edge correspondence for perfect matching faces
        void perfectPointMatch(const scalar absTol, const bool);

        //- Find point and edge correspondence for slaves being subdivision of
        //  master.
        void subDivisionMatch
        (
            const polyMesh& slaveMesh,
            const bool patchDivision,
            const scalar absTol
        );

public:

    //- Runtime type information
    ClassName("faceCoupleInfo");


    // Constructors

        //- Construct from two meshes and absolute tolerance.
        //  Finds out matches geometrically. No checking for nonsense match.
        //  Tolerance is absolute one so use with care.
        //  perfectMatch : each point/edge/face has corresponding point on other
        //                 side
        //                 if this is false then assumes slave is subdivision.
        //                 Matching then will work only for non-warped faces
        //                 since does nearest-to-face comparison with absTol.
        faceCoupleInfo
        (
            const polyMesh& mesh0,
            const polyMesh& mesh1,
            const scalar absTol,
            const bool perfectMatch
        );

        //- Construct from meshes and subset of mesh faces
        //  (i.e. indirectPrimitivePatch addressing)
        //  All faces in patch are considered matched (but don't have to be
        //  ordered)
        //  perfectMatch : each point/edge/face has corresponding point on other
        //                 side
        //  orderedFaces : faces in patch are ordered (so masterAddressing[i]
        //  matches slaveAddressing[i])
        //  patchDivision: faces in slave mesh that originate from the
        //  same master face have the same patch. Used by some triangulation
        //  methods.
        faceCoupleInfo
        (
            const polyMesh& masterMesh,
            const labelList& masterAddressing,
            const polyMesh& slaveMesh,
            const labelList& slaveAddressing,
            const scalar absTol,
            const bool perfectMatch,
            const bool orderedFaces,
            const bool patchDivision
        );


    //- Destructor
    ~faceCoupleInfo();



    // Member Functions

        //- Utility functions

            //- Get patch face labels
            static labelList faceLabels(const polyPatch&);

            //- Create Map from List
            static Map<label> makeMap(const labelList&);
            static Map<labelList> makeMap(const labelListList&);


        // Access

            //- Addressing engine for coupled faces on mesh0
            const indirectPrimitivePatch& masterPatch() const
            {
                return masterPatchPtr_();
            }

            //- Addressing engine for coupled faces on mesh1
            const indirectPrimitivePatch& slavePatch() const
            {
                return slavePatchPtr_();
            }

            //- Addressing engine for combined set of faces.
            const primitiveFacePatch& cutFaces() const
            {
                return cutFacesPtr_();
            }

            //- Points for combined set of faces.
            const pointField& cutPoints() const
            {
                return cutPoints_;
            }


        // Addressing from meshes to cut and vice versa.

            //- Master face for every face on cut. Will always be at least
            //  one but there might be multiple cut faces pointing to the same
            //  master
            const labelList& cutToMasterFaces() const
            {
                return cutToMasterFaces_;
            }
            const labelList& masterToCutPoints() const
            {
                return masterToCutPoints_;
            }

            const labelList& cutToSlaveFaces() const
            {
                return cutToSlaveFaces_;
            }
            const labelList& slaveToCutPoints() const
            {
                return slaveToCutPoints_;
            }

            //- From two cut points (original edge) to list of inserted
            //  points
            const edgeLookup& cutEdgeToPoints() const
            {
                return cutEdgeToPoints_;
            }

 };


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "faceCoupleInfoTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
