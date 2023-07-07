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
    Foam::blockMeshCartesianConfiguration

Description
    From a set of input surface geometry files and a set of configuration
    parameters, writes out a blockMeshDict configuration file. The mesh
    consists of a single block, aligned with Cartesian axes.

SourceFiles
    blockMeshCartesianConfiguration.C

\*---------------------------------------------------------------------------*/

#ifndef blockMeshCartesianConfiguration_H
#define blockMeshCartesianConfiguration_H

#include "blockMeshConfigurationBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class blockMeshCartesianConfiguration Declaration
\*---------------------------------------------------------------------------*/

class blockMeshCartesianConfiguration
:
    public blockMeshConfigurationBase
{
    // Private Data

        //- Number of cells in background mesh block
        Vector<label> nCells_;

        //- Refinement factor used to scale nCells
        const label refineFactor_;

        //- Clear the default patch entries for the background mesh
        const bool clearBoundary_;


    // Private Member Functions

        //- Calculate the parameters for the blockMeshDict file
        void calcBlockMeshDict();

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
        blockMeshCartesianConfiguration
        (
            const fileName& name,
            const fileName& dir,
            const Time& time,
            const meshingSurfaceList& surfaces,
            const Vector<label>& nCells,
            const label refineFactor,
            const HashTable<Pair<word>>& patchOpts,
            const bool clearBoundary
        );

        //- Disallow default bitwise copy construction
        blockMeshCartesianConfiguration
        (
            const blockMeshCartesianConfiguration&
        ) = delete;


    //- Destructor
    ~blockMeshCartesianConfiguration();


    // Member Functions

        //- Write the blockMeshDict
        void write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const blockMeshCartesianConfiguration&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
