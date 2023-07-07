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
    Foam::snappyHexMeshConfiguration

Description
    From a set of input surface geometry files and some configuration
    parameters, writes out a snappyHexMeshDict configuration file.

SourceFiles
    snappyHexMeshConfigurationI.H
    snappyHexMeshConfiguration.C

\*---------------------------------------------------------------------------*/

#ifndef snappyHexMeshConfiguration_H
#define snappyHexMeshConfiguration_H

#include "caseFileConfiguration.H"
#include "meshingSurfaceList.H"
#include "Tuple3.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class snappyHexMeshConfiguration Declaration
\*---------------------------------------------------------------------------*/

class snappyHexMeshConfiguration
:
    public caseFileConfiguration
{
    // Private Typedefs

        //- The surface type
        typedef meshingSurface::surfaceType surfaceType;


    // Private Data

        //- List of surfaces used for meshing
        const meshingSurfaceList& surfaces_;

        //- Refinement level applied across the snappyHexMeshDict file
        const label refinementLevel_;

        //- Level of refinement on speficied surfaces
        const List<Tuple2<word, label>>& surfaceLevels_;

        //- Refinement regions specified by surface and level
        const List<Tuple2<word, label>>& refinementRegions_;

        //- Refinement boxes with level of refinement
        const List<Tuple3<vector, vector, label>>& refinementBoxes_;

        //- Refinement distances with level of refinement
        const List<Tuple3<word, scalar, label>>& refinementDists_;

        //- Using implicit feature capturing?
        const bool implicitFeatures_;

        //- Number of layers at wall patches, default 0
        const label layers_;

        //- insidePoint parameter
        point insidePoint_;

        //- nCellsBetweenLevels parameter
        const label nCellsBetweenLevels_;


    // Private Member Functions

        //- Write switches
        void writeSnappySwitches();

        //- Write surface entry for geometry sub-dictionary
        void writeGeometrySurface(const label surfID);

        //- Write searchableBox entry for geometry sub-dictionary
        void writeSearchableBox(const label i);

        //- Write geometry sub-dictionary
        void writeSnappyGeometry();

        //- Write features list in castellatedMeshControls
        void writeFeatures();

        //- Write refinementSurfaces level entry
        void writeRefinementSurfacesLevel(const label rl);

        //- Write refinementSurfaces level entry using refinementLevel_
        void writeRefinementSurfacesLevel();

        //- Write refinementSurfaces level entry for a surface with name
        void writeRefinementSurfacesLevel(const word& name);

        //- Write refinementSurfaces wall patchInfo entry
        void writePatchInfo(const word& type, const word& group="");

        //- Write refinementSurfaces
        void writeRefinementSurfaces();

        //- Write a refinement surface region
        void writeRefinementSurfacesRegion
        (
            const word regionName,
            const List<word>& regions
        );

        //- Write refinement surface region information
        void writeRefinementSurfacesRegions
        (
            const wordList& inletRegions,
            const wordList& outletRegions
        );

        //- Write a refinement region
        void writeRefinementRegion(const word& name, const label level);

        //- Write refinementRegions
        void writeRefinementRegions();

        //- Write castellatedMeshControls
        void writeCastellatedMeshControls();

        //- Write snapControls
        void writeSnapControls();

        //- Write addLayersControls
        void writeAddLayersControls();

        //- Write writeFlags
        void writeWriteFlags();

public:

    // Constructors

        //- Construct from components
        snappyHexMeshConfiguration
        (
            const fileName& name,
            const fileName& dir,
            const Time& time,
            const meshingSurfaceList& surfaces,
            const label refinementLevel,
            const List<Tuple2<word, label>>& surfaceLevels,
            const List<Tuple2<word, label>>& refinementRegions,
            const List<Tuple3<vector, vector, label>>& refinementBoxes,
            const List<Tuple3<word, scalar, label>>& refinementDists,
            const bool implicitFeatures,
            const label layers,
            const point& insidePoint,
            const label nCellsBetweenLevels
        );

        //- Disallow default bitwise copy construction
        snappyHexMeshConfiguration(const snappyHexMeshConfiguration&) = delete;


    //- Destructor
    ~snappyHexMeshConfiguration();


    // Member Functions

        // Write the snappyHexMeshDict
        void write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const snappyHexMeshConfiguration&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
