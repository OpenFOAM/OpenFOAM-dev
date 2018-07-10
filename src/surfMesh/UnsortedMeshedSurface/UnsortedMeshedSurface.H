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
    Foam::UnsortedMeshedSurface

Description
    A surface geometry mesh, in which the surface zone information is
    conveyed by the 'zoneId' associated with each face.

    This form of surface description is particularly useful for reading in
    surface meshes from third-party formats (eg, obj, stl, gts, etc.). It
    can also be particularly useful for situations in which the surface
    many be adjusted in an arbitrary manner without worrying about needed
    to adjust the zone information (eg, surface refinement).

See also
    The Foam::MeshedSurface - which is organized as a surface mesh, but
    with independent zone information.

SourceFiles
    UnsortedMeshedSurface.C

\*---------------------------------------------------------------------------*/

#ifndef UnsortedMeshedSurface_H
#define UnsortedMeshedSurface_H

#include "MeshedSurface.H"
#include "surfZoneIdentifierList.H"
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
class IFstream;

template<class Face> class MeshedSurface;
template<class Face> class MeshedSurfaceProxy;
template<class Face> class UnsortedMeshedSurface;

/*---------------------------------------------------------------------------*\
                   Class UnsortedMeshedSurface Declaration
\*---------------------------------------------------------------------------*/

template<class Face>
class UnsortedMeshedSurface
:
    public MeshedSurface<Face>
{
    // friends - despite different face representationsx
    template<class Face2> friend class MeshedSurface;
    template<class Face2> friend class UnsortedMeshedSurface;
    friend class surfMesh;

private:

    // Private typedefs for convenience

        typedef MeshedSurface<Face>       ParentType;
        typedef MeshedSurface<Face>       FriendType;
        typedef MeshedSurfaceProxy<Face>  ProxyType;


    // Private Member Data

        //- The zone Id associated with each face
        labelList zoneIds_;

        //- Zone information (face ordering nFaces/startFace only used
        //  during reading and writing)
        List<surfZoneIdentifier> zoneToc_;


    // Private Member Functions

        //- Disable resize with value
        void resize(const label, const Face&);

        //- Disable setSize with value
        void setSize(const label, const Face&);


protected:

    // Protected Member functions

        //- Return non-const access to the zone Ids
        List<label>& storedZoneIds()
        {
            return zoneIds_;
        }

        //- Return non-const access to the zone table-of-contents
        List<surfZoneIdentifier>& storedZoneToc()
        {
            return zoneToc_;
        }

        //- Set new zones from faceMap
        virtual void remapFaces(const labelUList& faceMap);


public:

    // Public typedefs

        //- Face type used
        typedef Face FaceType;

        //- Runtime type information
        TypeName("UnsortedMeshedSurface");


    // Static

        //- Can we read this file format?
        static bool canReadType(const word& ext, const bool verbose=false);

        //- Can we read this file format?
        static bool canRead(const fileName&, const bool verbose=false);

        //- Can we write this file format?
        static bool canWriteType(const word& ext, const bool verbose=false);

        static wordHashSet readTypes();
        static wordHashSet writeTypes();


    // Constructors

        //- Construct null
        UnsortedMeshedSurface();

        //- Construct by transferring components
        //  (points, faces, zone ids, zone info).
        UnsortedMeshedSurface
        (
            const Xfer<pointField>&,
            const Xfer<List<Face>>&,
            const Xfer<List<label>>& zoneIds,
            const Xfer<surfZoneIdentifierList>&
        );

        //- Construct by transferring points, faces.
        //  Use zone information, or set single default zone
        UnsortedMeshedSurface
        (
            const Xfer<pointField>&,
            const Xfer<List<Face>>&,
            const labelUList& zoneSizes = labelUList(),
            const UList<word>& zoneNames = UList<word>()
        );

        //- Construct as copy
        UnsortedMeshedSurface(const UnsortedMeshedSurface<Face>&);

        //- Construct from a meshedSurface
        UnsortedMeshedSurface(const MeshedSurface<Face>&);

        //- Construct by transferring the contents from a UnsortedMeshedSurface
        UnsortedMeshedSurface(const Xfer<UnsortedMeshedSurface<Face>>&);

        //- Construct by transferring the contents from a meshedSurface
        UnsortedMeshedSurface(const Xfer<MeshedSurface<Face>>&);

        //- Construct from file name (uses extension to determine type)
        UnsortedMeshedSurface(const fileName&);

        //- Construct from file name (uses extension to determine type)
        UnsortedMeshedSurface(const fileName&, const word&);

        //- Construct from objectRegistry and a named surface
        UnsortedMeshedSurface(const Time&, const word& surfName="");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            UnsortedMeshedSurface,
            fileExtension,
            (
                const fileName& name
            ),
            (name)
        );


    // Selectors

        //- Select constructed from filename (explicit extension)
        static autoPtr<UnsortedMeshedSurface> New
        (
            const fileName&,
            const word& ext
        );

        //- Select constructed from filename (implicit extension)
        static autoPtr<UnsortedMeshedSurface> New(const fileName&);


    //- Destructor
    virtual ~UnsortedMeshedSurface();


    // Member Function Selectors

        declareMemberFunctionSelectionTable
        (
            void,
            UnsortedMeshedSurface,
            write,
            fileExtension,
            (
                const fileName& name,
                const UnsortedMeshedSurface<Face>& surf
            ),
            (name, surf)
        );

        //- Write to file
        static void write(const fileName&, const UnsortedMeshedSurface<Face>&);


    // Member Functions

        // Access

            //- The surface size is the number of faces
            label size() const
            {
                return ParentType::size();
            }

            //- Reset size of face and zone list
            void setSize(const label);

            //- Return const access to the zone ids
            const List<label>& zoneIds() const
            {
                return zoneIds_;
            }

            //- Return const access to the zone table-of-contents
            const List<surfZoneIdentifier>& zoneToc() const
            {
                return zoneToc_;
            }

            //- Sort faces according to zoneIds
            //  Returns a surfZoneList and sets faceMap to index within faces()
            //  (i.e. map from original,unsorted to sorted)
            surfZoneList sortedZones(labelList& faceMap) const;

            //- Set zones to 0 and set a single zone
            void setOneZone();

            //- Set zone ids and zones
            void setZones(const surfZoneList&);

            //- Set zone ids and zones
            void setZones(const labelUList& sizes, const UList<word>& names);

            //- Set zone ids and zones with default names
            void setZones(const labelUList& sizes);


        // Edit

            //- Clear all storage
            virtual void clear();

            //- Return new surface.
            //  Returns return pointMap, faceMap from subsetMeshMap
            UnsortedMeshedSurface subsetMesh
            (
                const labelHashSet& include,
                labelList& pointMap,
                labelList& faceMap
            ) const;

            //- Return new surface.
            UnsortedMeshedSurface subsetMesh
            (
                const labelHashSet& include
            ) const;

            //- Inherit reset from MeshedSurface<Face>
            using MeshedSurface<Face>::reset;

            //- Transfer components (points, faces, zone ids).
            virtual void reset
            (
                const Xfer<pointField>&,
                const Xfer<List<Face>>&,
                const Xfer<List<label>>& zoneIds
            );

            //- Transfer components (points, faces, zone ids).
            virtual void reset
            (
                const Xfer<List<point>>&,
                const Xfer<List<Face>>&,
                const Xfer<List<label>>& zoneIds
            );

            //- Transfer the contents of the argument and annul the argument
            void transfer(UnsortedMeshedSurface<Face>&);

            //- Transfer the contents of the argument and annul the argument
            void transfer(MeshedSurface<Face>&);

            //- Transfer contents to the Xfer container
            Xfer<UnsortedMeshedSurface<Face>> xfer();


        // Read

            //- Read from file. Chooses reader based on explicit extension
            bool read(const fileName&, const word& ext);

            //- Read from file. Chooses reader based on detected extension
            virtual bool read(const fileName&);


        // Write

            //- Generic write routine. Chooses writer based on extension.
            virtual void write(const fileName& name) const
            {
                write(name, *this);
            }

            //- Write to database
            void write(const Time&, const word& surfName="") const;


        // Member operators

            void operator=(const UnsortedMeshedSurface<Face>&);

            //- Conversion operator to MeshedSurfaceProxy
            operator MeshedSurfaceProxy<Face>() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "UnsortedMeshedSurface.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
