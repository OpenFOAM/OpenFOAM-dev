/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::boundaryCutter

Description
    Does modifications to boundary faces.

    Does
    - move boundary points
    - split boundary edges (multiple per edge if necessary)
    - face-centre decomposes boundary faces
    - diagonal split of boundary faces


SourceFiles
    boundaryCutter.C

\*---------------------------------------------------------------------------*/

#ifndef boundaryCutter_H
#define boundaryCutter_H

#include "Map.H"
#include "labelList.H"
#include "edge.H"
#include "typeInfo.H"
#include "labelPair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyTopoChange;
class mapPolyMesh;
class polyMesh;
class face;

/*---------------------------------------------------------------------------*\
                        Class boundaryCutter Declaration
\*---------------------------------------------------------------------------*/

class boundaryCutter
{
    // Private data

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Per edge sorted (start to end) list of points added.
        HashTable<labelList, edge, Hash<edge> > edgeAddedPoints_;

        //- Per face the mid point added.
        Map<label> faceAddedPoint_;


    // Private Member Functions

        //- Get patch and zone info for face
        void getFaceInfo
        (
            const label faceI,
            label& patchID,
            label& zoneID,
            label& zoneFlip
        ) const;

        //- Add cuts of edges to face
        face addEdgeCutsToFace(const label faceI, const Map<labelList>&) const;

        //- Splits faces with multiple cut edges. Return true if anything split.
        bool splitFace
        (
            const label faceI,
            const Map<point>& pointToPos,
            const Map<labelList>& edgeToAddedPoints,
            polyTopoChange& meshMod
        ) const;

        //- Add/modify faceI for new vertices.
        void addFace
        (
            const label faceI,
            const face& newFace,

            bool& modifiedFace,     // have we already 'used' faceI?
            polyTopoChange& meshMod
        ) const;


        //- Disallow default bitwise copy construct
        boundaryCutter(const boundaryCutter&);

        //- Disallow default bitwise assignment
        void operator=(const boundaryCutter&);

public:

    //- Runtime type information
    ClassName("boundaryCutter");

    // Constructors

        //- Construct from mesh
        boundaryCutter(const polyMesh& mesh);


    //- Destructor
    ~boundaryCutter();


    // Member Functions

        // Edit

            //- Do actual cutting with cut description. Inserts mesh changes
            //  into meshMod.
            //  pointToPos : new position for selected points
            //  edgeToCuts : per edge set of points that need to be introduced
            //  faceToSplit : per face the diagonal split
            //  faceToFeaturePoint : per face the feature point. Triangulation
            //                       around this feature point.
            void setRefinement
            (
                const Map<point>& pointToPos,
                const Map<List<point> >& edgeToCuts,
                const Map<labelPair>& faceToSplit,
                const Map<point>& faceToFeaturePoint,
                polyTopoChange& meshMod
            );

            //- Force recalculation of locally stored data on topological change
            void updateMesh(const mapPolyMesh&);


        // Access

            //- Per edge a sorted list (start to end) of added points.
            const HashTable<labelList, edge, Hash<edge> >& edgeAddedPoints()
             const
            {
                return edgeAddedPoints_;
            }

            //- Per face the mid point added.
            const Map<label>& faceAddedPoint() const
            {
                return faceAddedPoint_;
            }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
