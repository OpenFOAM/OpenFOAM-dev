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
    Foam::triSurface

Description
    Triangulated surface description with patch information.

SourceFiles
    triSurface.C

\*---------------------------------------------------------------------------*/

#ifndef triSurface_H
#define triSurface_H

#include "PrimitivePatch.H"
#include "pointField.H"
#include "labelledTri.H"
#include "boolList.H"
#include "geometricSurfacePatchList.H"
#include "surfacePatchList.H"
#include "triFaceList.H"
#include "triadField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Time;
class IFstream;


// Forward declaration of friend functions and operators

class triSurface;

Ostream& operator<<(Ostream&, const triSurface&);


/*---------------------------------------------------------------------------*\
                           Class triSurface Declaration
\*---------------------------------------------------------------------------*/

class triSurface
:
    public PrimitivePatch<labelledTri, ::Foam::List, pointField, point>
{
    // Private typedefs

    //- Typedefs for convenience
        typedef labelledTri Face;
        typedef PrimitivePatch
        <
            labelledTri,
            ::Foam::List,
            pointField,
            point
        >
        ParentType;


    // Private data

        //- The number of bytes in the STL header
        static const int STLheaderSize = 80;

        //- Patch information (face ordering nFaces/startFace only used
        //  during reading and writing)
        geometricSurfacePatchList patches_;


    // Demand driven private data.

        //- Edge-face addressing (sorted)
        mutable labelListList* sortedEdgeFacesPtr_;

        //- Label of face that 'owns' edge (i.e. e.vec() is righthanded walk
        //  along face)
        mutable labelList* edgeOwnerPtr_;


    // Private Member Functions

        //- Calculate sorted edgeFaces
        void calcSortedEdgeFaces() const;

        //- Calculate owner
        void calcEdgeOwner() const;

        //- Sort faces according to region. Returns patch list
        //  and sets faceMap to index of labelledTri inside *this.
        surfacePatchList calcPatches(labelList& faceMap) const;

        //- Sets default values for patches
        void setDefaultPatches();

        //- Function to stitch the triangles by removing duplicate points.
        //  Returns true if any points merged
        bool stitchTriangles
        (
            const scalar tol = small,
            const bool verbose = false
        );

        scalar pointNormalWeight
        (
            const triFace& f,
            const label pi,
            const vector& fa,
            const pointField& points
        ) const;

        //- Return the surface point normals
        tmp<vectorField> weightedPointNormals() const;

        //- Return the curvature of surface at the points
        tmp<triadField> pointCoordSys(const vectorField& pointNormals) const;


        //- Read in Foam format
        bool read(Istream&);

        //- Generic read routine. Chooses reader based on extension.
        bool read(const fileName&, const word& ext, const bool check = true);

        bool readSTL(const fileName&);
        bool readSTLASCII(const fileName&);
        bool readSTLBINARY(const fileName&);
        bool readGTS(const fileName&);
        bool readOBJ(const fileName&);
        bool readOFF(const fileName&);
        bool readTRI(const fileName&);
        bool readAC(const fileName&);
        bool readNAS(const fileName&);
        bool readVTK(const fileName&);

        //- Generic write routine. Chooses writer based on extension.
        void write(const fileName&, const word& ext, const bool sort) const;

        //- Write to Ostream in ASCII STL format.
        //  Each region becomes 'solid' 'endsolid' block.
        void writeSTLASCII(const bool writeSorted, Ostream&) const;

        //- Write to std::ostream in BINARY STL format
        void writeSTLBINARY(std::ostream&) const;

        //- Write to Ostream in GTS (Gnu Tri Surface library)
        //  format.
        void writeGTS(const bool writeSorted, Ostream&) const;

        //- Write to Ostream in OBJ (Lightwave) format.
        //  writeSorted=true: sort faces acc. to region and write as single
        //  group. =false: write in normal order.
        void writeOBJ(const bool writeSorted, Ostream&) const;

        //- Write to Ostream in OFF (Geomview) format.
        //  writeSorted=true: sort faces acc. to region and write as single
        //  group. =false: write in normal order.
        void writeOFF(const bool writeSorted, Ostream&) const;

        //- Write to VTK legacy format.
        void writeVTK(const bool writeSorted, Ostream&) const;

        //- Write to Ostream in TRI (AC3D) format
        //  Ac3d .tri format (unmerged triangle format)
        void writeTRI(const bool writeSorted, Ostream&) const;

        //- Write to Ostream in SMESH (tetgen) format
        void writeSMESH(const bool writeSorted, Ostream&) const;

        //- Write to Ostream in AC3D format. Always sorted by patch.
        void writeAC(Ostream&) const;


    // Static private functions

        //- Convert faces to labelledTri. All get same region.
        static List<labelledTri> convertToTri
        (
            const faceList&,
            const label defaultRegion = 0
        );

        //- Convert triFaces to labelledTri. All get same region.
        static List<labelledTri> convertToTri
        (
            const triFaceList&,
            const label defaultRegion = 0
        );

        //- Helper function to print triangle info
        static void printTriangle
        (
            Ostream&,
            const Foam::string& pre,
            const labelledTri&,
            const pointField&
        );

        //- Read non-comment line
        static string getLineNoComment(IFstream&);


protected:

    // Protected Member Functions

        //- Non-const access to global points
        pointField& storedPoints()
        {
            return const_cast<pointField&>(ParentType::points());
        }

        //- Non-const access to the faces
        List<Face>& storedFaces()
        {
            return static_cast<List<Face>&>(*this);
        }


public:

    // Public typedefs

        //- Placeholder only, but do not remove - it is needed for GeoMesh
        typedef bool BoundaryMesh;


        //- Runtime type information
        ClassName("triSurface");


    // Static

        //- Name of triSurface directory to use.
        static fileName triSurfInstance(const Time&);


    // Constructors

        //- Construct null
        triSurface();

        //- Construct from triangles, patches, points.
        triSurface
        (
            const List<labelledTri>&,
            const geometricSurfacePatchList&,
            const pointField&
        );

        //- Construct from triangles, patches, points. Reuse storage.
        triSurface
        (
            List<labelledTri>&,
            const geometricSurfacePatchList&,
            pointField&,
            const bool reuse
        );

        //- Construct from triangles, patches, points.
        triSurface
        (
            const Xfer<List<labelledTri>>&,
            const geometricSurfacePatchList&,
            const Xfer<List<point>>&
        );

        //- Construct from triangles, points. Set patchnames to default.
        triSurface(const List<labelledTri>&, const pointField&);

        //- Construct from triangles, points. Set region to 0 and default
        //  patchName.
        triSurface(const triFaceList&, const pointField&);

        //- Construct from file name (uses extension to determine type)
        triSurface(const fileName&);

        //- Construct from Istream
        triSurface(Istream&);

        //- Construct from objectRegistry
        triSurface(const Time& d);

        //- Construct as copy
        triSurface(const triSurface&);


    //- Destructor
    ~triSurface();

        void clearOut();

        void clearTopology();

        void clearPatchMeshAddr();


    // Member Functions

        // Access

            const geometricSurfacePatchList& patches() const
            {
                return patches_;
            }

            geometricSurfacePatchList& patches()
            {
                return patches_;
            }

            //- Return edge-face addressing sorted (for edges with more than
            //  2 faces) according to the angle around the edge.
            //  Orientation is anticlockwise looking from
            //  edge.vec(localPoints())
            const labelListList& sortedEdgeFaces() const;

            //- If 2 face neighbours: label of face where ordering of edge
            //  is consistent with righthand walk.
            //  If 1 neighbour: label of only face.
            //  If >2 neighbours: undetermined.
            const labelList& edgeOwner() const;


        // Edit

            //- Move points
            virtual void movePoints(const pointField&);

            //- Scale points. A non-positive factor is ignored
            virtual void scalePoints(const scalar);

            //- Check/remove duplicate/degenerate triangles
            void checkTriangles(const bool verbose);

            //- Check triply (or more) connected edges.
            void checkEdges(const bool verbose);

            //- Remove non-valid triangles
            void cleanup(const bool verbose);

            //- Fill faceZone with currentZone for every face reachable
            //  from facei without crossing edge marked in borderEdge.
            //  Note: faceZone has to be sized nFaces before calling this fun.
            void markZone
            (
                const boolList& borderEdge,
                const label facei,
                const label currentZone,
                labelList& faceZone
            ) const;

            //- (size and) fills faceZone with zone of face. Zone is area
            //  reachable by edge crossing without crossing borderEdge
            //  (bool for every edge in surface). Returns number of zones.
            label markZones
            (
                const boolList& borderEdge,
                labelList& faceZone
            ) const;

            //- 'Create' sub mesh, including only faces for which
            //  boolList entry is true
            //  Sets: pointMap: from new to old localPoints
            //        faceMap: new to old faces
            void subsetMeshMap
            (
                const boolList& include,
                labelList& pointMap,
                labelList& faceMap
            ) const;

            //- Return new surface. Returns pointMap, faceMap from
            //  subsetMeshMap
            triSurface subsetMesh
            (
                const boolList& include,
                labelList& pointMap,
                labelList& faceMap
            ) const;


        // Conversion

            //- Return the list of triangles as a faceList
            faceList faces() const;


        // Analysis

            //- Return the curvature of surface at the points
            tmp<scalarField> curvature() const;


        // Write

            //- Write to Ostream in simple FOAM format
            void write(Ostream&) const;

            //- Generic write routine. Chooses writer based on extension.
            void write(const fileName&, const bool sortByRegion = false) const;

            //- Write to database
            void write(const Time&) const;

            //- Write some statistics
            void writeStats(Ostream&) const;


    // Member operators

        void operator=(const triSurface&);


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const triSurface&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
