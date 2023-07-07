/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::meshingSurface

Description
    Attributes of a surface geometry file (e.g. OBJ, STL) that are used in
    the configuration of mesh input files, (e.g. blockMeshDict,
    snappyHexMeshDict).

SourceFiles
    meshingSurface.C

\*---------------------------------------------------------------------------*/

#ifndef meshingSurface_H
#define meshingSurface_H

#include "Time.H"
#include "boundBox.H"
#include "NamedEnum.H"
#include "triSurfaceMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class meshingSurface Declaration
\*---------------------------------------------------------------------------*/

class meshingSurface
{
public:

    // Static Data Members

        //- Type of surface
        enum class surfaceType
        {
            wall,         // Solid wall boundary
            external,     // External domain boundary
            cellZone,     // Surface defining a cellZone
            rotatingZone, // Surface defining a rotatingZone
            baffle        // Zero-thickness wall
        };


private:

    // Private Data

        //- Path to surface file
        fileName path_;

        //- Surface file name
        fileName file_;

        //- Surfaces names, without file extensions
        word name_;

        //- Type of surface: wall, external, cellZone, rotatingZone, baffle
        surfaceType type_;

        //- Bounding boxes of the surface
        boundBox boundBox_;

        //- Is the surface closed?
        bool closed_;

        //- Number of disconnected parts
        label nParts_;

        //- List of regions
        List<word> regions_;

        //- List of inlet regions
        List<word> inletRegions_;

        //- List of outlet regions
        List<word> outletRegions_;


    // Private Member Functions

        //- Returns the number of disconnected regions in a surface
        label nSurfaceParts(const triSurfaceMesh& surf);


public:

    // Static Data Members

        //- Enumeration names for surfaceTypes
        static const NamedEnum<surfaceType, 5> surfaceTypeNames;


    // Constructors

        //- Construct null
        meshingSurface();

        //- Construct from components
        meshingSurface(const fileName& file, const Time& time);

        //- Disallow default bitwise copy construction
        meshingSurface(const meshingSurface&) = delete;


    //- Destructor
    ~meshingSurface();


    // Member Functions

        //- Path to the surface file
        const fileName& path() const
        {
            return path_;
        }

        //- Surface file name
        const fileName& file() const
        {
            return file_;
        }

        //- Surface name, without file extension
        const word& name() const
        {
            return name_;
        }

        //- Surface type
        const surfaceType& type() const
        {
            return type_;
        }

        //- Return non-const access to the surface type
        surfaceType& type()
        {
            return type_;
        }

        //- Surface bounding box
        const boundBox& bb() const
        {
            return boundBox_;
        }

        //- Is the surface closed?
        bool closed() const
        {
            return closed_;
        }

        //- Return the number of disconnected surfaces
        label nParts() const
        {
            return nParts_;
        }

        //- Surface geometry regions
        const List<word>& regions() const
        {
            return regions_;
        }

        //- Inlet regions
        const List<word>& inletRegions() const
        {
            return inletRegions_;
        }

        //- Inlet regions
        List<word>& inletRegions()
        {
            return inletRegions_;
        }

        //- Outlet regions
        const List<word>& outletRegions() const
        {
            return outletRegions_;
        }

        //- Outlet regions
        List<word>& outletRegions()
        {
            return outletRegions_;
        }

        // Check the file extension of surface geometry files
        static bool isSurfaceExt(const fileName& file);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const meshingSurface&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
