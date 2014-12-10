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
    Foam::cellSplitter

Description
    Does pyramidal decomposition of selected cells. So all faces will become
    base of pyramid with as top a user-supplied point (usually the cell centre)

SourceFiles
    cellSplitter.C

\*---------------------------------------------------------------------------*/

#ifndef cellSplitter_H
#define cellSplitter_H

#include "Map.H"
#include "edge.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyTopoChange;
class mapPolyMesh;
class polyMesh;


/*---------------------------------------------------------------------------*\
                           Class cellSplitter Declaration
\*---------------------------------------------------------------------------*/

class cellSplitter
{
    // Private data

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Per cell the mid point added.
        Map<label> addedPoints_;


    // Private Member Functions

        //- Get patch and zone info for face
        void getFaceInfo
        (
            const label faceI,
            label& patchID,
            label& zoneID,
            label& zoneFlip
        ) const;

        //- Find the new owner (if any) of the face.
        label newOwner
        (
            const label faceI,
            const Map<labelList>& cellToCells
        ) const;

        //- Find the new neighbour (if any) of the face.
        label newNeighbour
        (
            const label faceI,
            const Map<labelList>& cellToCells
        ) const;

        //- Disallow default bitwise copy construct
        cellSplitter(const cellSplitter&);

        //- Disallow default bitwise assignment
        void operator=(const cellSplitter&);

public:

    //- Runtime type information
    ClassName("cellSplitter");

    // Constructors

        //- Construct from mesh
        cellSplitter(const polyMesh& mesh);


    //- Destructor
    ~cellSplitter();


    // Member Functions

        // Edit

            //- Insert mesh changes into meshMod.
            //  cellToMidPoint : cell to cut and position of its new midpoint
            void setRefinement
            (
                const Map<point>& cellToMidPoint,
                polyTopoChange& meshMod
            );

            //- Force recalculation of locally stored data on topological change
            void updateMesh(const mapPolyMesh&);


        // Access

            //- Per cell the mid point added.
            const Map<label>& addedPoints() const
            {
                return addedPoints_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
