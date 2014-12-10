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
    Foam::MeshedSurfaceProxy

Description
    A proxy for writing MeshedSurface, UnsortedMeshedSurface and surfMesh
    to various file formats.

SourceFiles
    MeshedSurfaceProxy.C
    MeshedSurfaceProxyCore.C

\*---------------------------------------------------------------------------*/

#ifndef MeshedSurfaceProxy_H
#define MeshedSurfaceProxy_H

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

template<class Face> class MeshedSurface;

/*---------------------------------------------------------------------------*\
                     Class MeshedSurfaceProxy Declaration
\*---------------------------------------------------------------------------*/

template<class Face>
class MeshedSurfaceProxy
:
    public fileFormats::surfaceFormatsCore
{
    //- Private data

        const pointField& points_;

        const List<Face>& faces_;

        const List<surfZone>& zones_;

        const List<label>& faceMap_;


public:

    // Public typedefs

        //- Face type used
        typedef Face FaceType;


    // Static

        //- Runtime type information
        ClassName("MeshedSurfaceProxy");

        //- The file format types that can be written via MeshedSurfaceProxy
        static wordHashSet writeTypes();

        //- Can this file format type be written via MeshedSurfaceProxy?
        static bool canWriteType(const word& ext, const bool verbose=false);


    // Constructors

        //- Construct from component references
        MeshedSurfaceProxy
        (
            const pointField&,
            const List<Face>&,
            const List<surfZone>& = List<surfZone>(),
            const List<label>& faceMap = List<label>()
        );


    //- Destructor
    virtual ~MeshedSurfaceProxy();


    // Member Function Selectors

        declareMemberFunctionSelectionTable
        (
            void,
            MeshedSurfaceProxy,
            write,
            fileExtension,
            (
                const fileName& name,
                const MeshedSurfaceProxy<Face>& surf
            ),
            (name, surf)
        );

        //- Write to file
        static void write(const fileName&, const MeshedSurfaceProxy<Face>&);


    // Member Functions

        // Access

            //- Return const access to the points
            inline const pointField& points() const
            {
                return points_;
            }

            //- Return const access to the faces
            inline const List<Face>& faces() const
            {
                return faces_;
            }

            //- Const access to the surface zones.
            //  If zones are defined, they must be contiguous and cover the
            //  entire surface
            inline const List<surfZone>& surfZones() const
            {
                return zones_;
            }

            //- Const access to the faceMap, zero-sized when unused
            inline const List<label>& faceMap() const
            {
                return faceMap_;
            }

            //- Use faceMap?
            inline bool useFaceMap() const
            {
                return faceMap_.size() == faces_.size();
            }

        // Write

            //- Generic write routine. Chooses writer based on extension.
            virtual void write(const fileName& name) const
            {
                write(name, *this);
            }

            //- Write to database
            virtual void write(const Time&, const word& surfName = "") const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "MeshedSurfaceProxy.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
