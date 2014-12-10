/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::MeshedSurface

Description
    A surface geometry mesh with zone information, not to be confused with
    the similarly named surfaceMesh, which actually refers to the cell faces
    of a volume mesh.

    A MeshedSurface can have zero or more surface zones (roughly equivalent
    to faceZones for a polyMesh). If surface zones are defined, they must
    be contiguous and cover all of the faces.

    The MeshedSurface is intended for surfaces from a variety of sources.
    - A set of points and faces without any surface zone information.
    - A set of points and faces with randomly ordered zone information.
      This could arise, for example, from reading external file formats
      such as STL, etc.

SourceFiles
    MeshedSurface.C

\*---------------------------------------------------------------------------*/

#ifndef MeshedSurface_H
#define MeshedSurface_H

#include "PrimitivePatch.H"
#include "PatchTools.H"
#include "pointField.H"
#include "face.H"
#include "triFace.H"

#include "surfZoneList.H"
#include "surfaceFormatsCore.H"
#include "runTimeSelectionTables.H"
#include "memberFunctionSelectionTables.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class Time;
class surfMesh;
class polyBoundaryMesh;

template<class Face> class MeshedSurface;
template<class Face> class MeshedSurfaceProxy;
template<class Face> class UnsortedMeshedSurface;

/*---------------------------------------------------------------------------*\
                      Class MeshedSurface Declaration
\*---------------------------------------------------------------------------*/

template<class Face>
class MeshedSurface
:
    public PrimitivePatch<Face, ::Foam::List, pointField, point>,
    public fileFormats::surfaceFormatsCore
{
    // friends - despite different face representationsx
    template<class Face2> friend class MeshedSurface;
    template<class Face2> friend class UnsortedMeshedSurface;
    friend class surfMesh;


private:

    // Private typedefs for convenience

        typedef PrimitivePatch
        <
            Face,
            ::Foam::List,
            pointField,
            point
        >
        ParentType;

        typedef UnsortedMeshedSurface<Face>  FriendType;
        typedef MeshedSurfaceProxy<Face>     ProxyType;


    // Private Member Data

        //- Zone information
        // (face ordering nFaces/startFace only used during reading/writing)
        List<surfZone> zones_;


protected:

    // Protected Member functions

        //- Transfer points/zones and transcribe face -> triFace
        void transcribe(MeshedSurface<face>&);

        //- basic sanity check on zones
        void checkZones();

        //- Non-const access to global points
        pointField& storedPoints()
        {
            return const_cast<pointField&>(ParentType::points());
        }

        //- Non-const access to the faces
        List<Face>& storedFaces()
        {
            return static_cast<List<Face> &>(*this);
        }

        //- Non-const access to the zones
        surfZoneList& storedZones()
        {
            return zones_;
        }

        //- sort faces by zones and store sorted faces
        void sortFacesAndStore
        (
            const Xfer<List<Face> >& unsortedFaces,
            const Xfer<List<label> >& zoneIds,
            const bool sorted
        );

        //- Set new zones from faceMap
        virtual void remapFaces(const labelUList& faceMap);


public:

    // Public typedefs

        //- Face type used
        typedef Face FaceType;

        //- Runtime type information
        ClassName("MeshedSurface");


    // Static

        //- Face storage only handles triangulated faces
        inline static bool isTri();

        //- Can we read this file format?
        static bool canRead(const fileName&, const bool verbose=false);

        //- Can we read this file format?
        static bool canReadType(const word& ext, const bool verbose=false);

        //- Can we write this file format?
        static bool canWriteType(const word& ext, const bool verbose=false);

        static wordHashSet readTypes();
        static wordHashSet writeTypes();


    // Constructors

        //- Construct null
        MeshedSurface();

        //- Construct by transferring components (points, faces, zones).
        MeshedSurface
        (
            const Xfer<pointField>&,
            const Xfer<List<Face> >&,
            const Xfer<surfZoneList>&
        );

        //- Construct by transferring components (points, faces).
        //  Use zone information if available
        MeshedSurface
        (
            const Xfer<pointField>&,
            const Xfer<List<Face> >&,
            const labelUList& zoneSizes = labelUList(),
            const UList<word>& zoneNames = UList<word>()
        );

        //- Construct as copy
        MeshedSurface(const MeshedSurface&);

        //- Construct from a UnsortedMeshedSurface
        MeshedSurface(const UnsortedMeshedSurface<Face>&);

        //- Construct from a boundary mesh with local points/faces
        MeshedSurface
        (
            const polyBoundaryMesh&,
            const bool globalPoints=false
        );

        //- Construct from a surfMesh
        MeshedSurface(const surfMesh&);

        //- Construct by transferring the contents from a UnsortedMeshedSurface
        MeshedSurface(const Xfer<UnsortedMeshedSurface<Face> >&);

        //- Construct by transferring the contents from a MeshedSurface
        MeshedSurface(const Xfer<MeshedSurface<Face> >&);

        //- Construct from file name (uses extension to determine type)
        MeshedSurface(const fileName&);

        //- Construct from file name (uses extension to determine type)
        MeshedSurface(const fileName&, const word& ext);

        //- Construct from database
        MeshedSurface(const Time&, const word& surfName="");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            MeshedSurface,
            fileExtension,
            (
                const fileName& name
            ),
            (name)
        );


    // Selectors

        //- Select constructed from filename (explicit extension)
        static autoPtr<MeshedSurface> New
        (
            const fileName&,
            const word& ext
        );

        //- Select constructed from filename (implicit extension)
        static autoPtr<MeshedSurface> New(const fileName&);


    //- Destructor
    virtual ~MeshedSurface();


    // Member Function Selectors

        declareMemberFunctionSelectionTable
        (
            void,
            UnsortedMeshedSurface,
            write,
            fileExtension,
            (
                const fileName& name,
                const MeshedSurface<Face>& surf
            ),
            (name, surf)
        );

        //- Write to file
        static void write(const fileName&, const MeshedSurface<Face>&);


    // Member Functions

        // Access

            //- The surface size is the number of faces
            label size() const
            {
                return ParentType::size();
            }

            //- Return const access to the faces
            inline const List<Face>& faces() const
            {
                return static_cast<const List<Face> &>(*this);
            }

            //- Const access to the surface zones.
            //  If zones are defined, they must be contiguous and cover the
            //  entire surface
            const List<surfZone>& surfZones() const
            {
                return zones_;
            }

            //- Add surface zones
            virtual void addZones
            (
                const UList<surfZone>&,
                const bool cullEmpty=false
            );

            //- Add surface zones
            virtual void addZones
            (
                const labelUList& sizes,
                const UList<word>& names,
                const bool cullEmpty=false
            );

            //- Add surface zones
            virtual void addZones
            (
                const labelUList& sizes,
                const bool cullEmpty=false
            );

            //- Remove surface zones
            virtual void removeZones();


        // Edit

            //- Clear all storage
            virtual void clear();

            //- Move points
            virtual void movePoints(const pointField&);

            //- Scale points. A non-positive factor is ignored
            virtual void scalePoints(const scalar);

            //- Reset primitive data (points, faces and zones)
            //  Note, optimized to avoid overwriting data (with Xfer::null)
            virtual void reset
            (
                const Xfer<pointField >& points,
                const Xfer<List<Face> >& faces,
                const Xfer<surfZoneList>& zones
            );

            //- Reset primitive data (points, faces and zones)
            //  Note, optimized to avoid overwriting data (with Xfer::null)
            virtual void reset
            (
                const Xfer<List<point> >& points,
                const Xfer<List<Face> >& faces,
                const Xfer<surfZoneList >& zones
            );

            //- Remove invalid faces
            virtual void cleanup(const bool verbose);

            virtual bool stitchFaces
            (
                const scalar tol=SMALL,
                const bool verbose=false
            );

            virtual bool checkFaces
            (
                const bool verbose=false
            );

            //- Triangulate in-place, returning the number of triangles added
            virtual label triangulate();

            //- Triangulate in-place, returning the number of triangles added
            //  and setting a map of original face Ids.
            //  The faceMap is zero-sized when no triangulation was done.
            virtual label triangulate(List<label>& faceMap);


            //- Return new surface.
            //  Returns return pointMap, faceMap from subsetMeshMap
            MeshedSurface subsetMesh
            (
                const labelHashSet& include,
                labelList& pointMap,
                labelList& faceMap
            ) const;

            //- Return new surface.
            MeshedSurface subsetMesh
            (
                const labelHashSet& include
            ) const;

            //- Transfer the contents of the argument and annul the argument
            void transfer(MeshedSurface<Face>&);

            //- Transfer the contents of the argument and annul the argument
            void transfer(UnsortedMeshedSurface<Face>&);

            //- Transfer contents to the Xfer container
            Xfer<MeshedSurface<Face> > xfer();


        // Read

            //- Read from file. Chooses reader based on explicit extension
            bool read(const fileName&, const word& ext);

            //- Read from file. Chooses reader based on detected extension
            virtual bool read(const fileName&);


        // Write

            void writeStats(Ostream& os) const;

            //- Generic write routine. Chooses writer based on extension.
            virtual void write(const fileName& name) const
            {
                write(name, *this);
            }

            //- Write to database
            void write(const Time&, const word& surfName="") const;


    // Member operators

        void operator=(const MeshedSurface<Face>&);

        //- Conversion operator to MeshedSurfaceProxy
        operator MeshedSurfaceProxy<Face>() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Specialization for holding triangulated information
template<>
inline bool MeshedSurface<triFace>::isTri()
{
    return true;
}


//- Specialization for holding triangulated information
template<>
inline label MeshedSurface<triFace>::triangulate()
{
    return 0;
}


//- Specialization for holding triangulated information
template<>
inline label MeshedSurface<triFace>::triangulate(List<label>& faceMap)
{
    if (&faceMap)
    {
        faceMap.clear();
    }

    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "MeshedSurface.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
