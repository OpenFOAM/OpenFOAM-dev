/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2024 OpenFOAM Foundation
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
    Foam::blockMeshCylindricalConfiguration

Description
    From a set of input surface geometry files and a set of configuration
    parameters, writes out a blockMeshDict configuration file. The mesh
    consists of a single block, aligned with cylindrical coordinates about
    the z-axis.

SourceFiles
    blockMeshCylindricalConfiguration.C

\*---------------------------------------------------------------------------*/

#ifndef blockMeshCylindricalConfiguration_H
#define blockMeshCylindricalConfiguration_H

#include "blockMeshConfigurationBase.H"
#include "meshingSurfaceList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class blockMeshCylindricalConfiguration Declaration
\*---------------------------------------------------------------------------*/

class blockMeshCylindricalConfiguration
:
    public blockMeshConfigurationBase
{
    // Private Data

        //- Bounding box for the rotatingZone surfaces
        boundBox rzbb_;

        //- Half width of the core box
        scalar radBox_;

        //- Number of cells in background mesh = (boxCells radialCells zCells)
        Vector<label> nCells_;

        //- Refinement factor used to scale nCells
        label refineFactor_;

        //- Clear the default patch entries for the background mesh
        bool clearBoundary_;


    // Private Member Functions

        //- Inflate a bounding box by a scaling vector
        void bbInflate(boundBox& bb, const vector& s);

        //- Are the geometry bounds aligned with the z-axis
        bool isBoundBoxOnZaxis();

        //- Calculate the parameters for the blockMeshDict file
        void calcBlockMeshDict
        (
            const bool& boundsOpt,
            const bool& rotatingZonesOpt
        );

        //- Write backgroundMesh sub-dictionary
        void writeBackgroundMesh();

        //- Write the defaultPatch entry
        void writeDefaultPatch();

        //- Write a patch in the boundary sub-dictionary
        void writePatch
        (
            const word& name,
            const word& type,
            const string& face
        );

        //- Write the boundary sub-dictionary
        void writeBoundary();

        //- Write the geometry sub-dictionary
        void writeGeometry();

        //- Write a projected vertex entry in the vertices list
        void writeProjectedVertex
        (
            const word& x,
            const word& y,
            const word& z,
            const word& surface
        );

        //- Write vertices list
        void writeVertices();

        //- Write blocks sub-dictionary
        void writeBlocks();

        //- Write edges list
        void writeEdges();

        //- Write mergePatchPairs
        void writeMergePatchPairs();


public:

    // Static Data Members

        //- Default patch names for the background mesh
        static const List<word> patches;


    // Constructors

        //- Construct from components
        blockMeshCylindricalConfiguration
        (
            const fileName& name,
            const fileName& dir,
            const Time& time,
            const meshingSurfaceList& surfaces,
            const bool& boundsOpt,
            const Vector<label>& nCells,
            const label refineFactor,
            const HashTable<Pair<word>>& patchOpts,
            const bool clearBoundary
        );

        //- Disallow default bitwise copy construction
        blockMeshCylindricalConfiguration
        (
            const blockMeshCylindricalConfiguration&
        ) = delete;


    //- Destructor
    ~blockMeshCylindricalConfiguration();


    // Member Functions

        // Write the blockMeshDict
        void write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const blockMeshCylindricalConfiguration&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
