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
    Foam::removePoints

Description
    Removes selected points from mesh and updates faces using these
    points.

SourceFiles
    removePoints.C

\*---------------------------------------------------------------------------*/

#ifndef removePoints_H
#define removePoints_H

#include "typeInfo.H"
#include "boolList.H"
#include "pointField.H"
#include "faceList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class polyTopoChange;
class mapPolyMesh;
class face;

/*---------------------------------------------------------------------------*\
                           Class removePoints Declaration
\*---------------------------------------------------------------------------*/

class removePoints
{

    // Private data

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Whether undoable
        const bool undoable_;

        //- If undoable: deleted points
        pointField savedPoints_;

        //- If undoable: per stored face the original mesh face label
        labelList savedFaceLabels_;

        //- If undoable: per stored face the vertices. Negative indices
        //  refer to deletedPoints_
        faceList savedFaces_;


    // Private Member Functions

        //- Change the vertices of the face whilst keeping everything else
        //  (patch, zone) the same.
        void modifyFace
        (
            const label facei,
            const face&,
            polyTopoChange&
        ) const;


        //- Disallow default bitwise copy construct
        removePoints(const removePoints&);

        //- Disallow default bitwise assignment
        void operator=(const removePoints&);

public:

    //- Runtime type information
    ClassName("removePoints");


    // Constructors

        //- Construct from mesh
        removePoints(const polyMesh& mesh, const bool undoable = false);


    // Member Functions

            //- If undoable: affected face labels. Already restored faces
            //  will be -1.
            const labelList& savedFaceLabels() const
            {
                return savedFaceLabels_;
            }


        // Helper functions

            //- Mark in pointCanBeDeleted the points that can be deleted
            //  (parallel synchronised) and returns the global number of these
            //  points. (this number is the global number before synchronisation
            //  so might be off!)
            //  A point can be deleted if
            //  - it is not used by any edge.
            //  or
            //  - is used by only two edges
            //  - these two edges are sufficiently in line (cos > minCos)
            //  - all processors agree that point can be deleted.
            label countPointUsage
            (
                const scalar minCos,
                boolList& pointCanBeDeleted
            ) const;

        // Topology changes

            //- Play commands into polyTopoChange to remove points. Gets
            //  boolList (output of countPointUsage) as input.
            //  Does no check for whether resulting face is legal.
            //  Since pointCanBeDeleted is synced all coupled faces should
            //  decide the same.
            void setRefinement(const boolList&, polyTopoChange&);

            //- Force recalculation of locally stored data on topological change
            void updateMesh(const mapPolyMesh&);

            //- Given set of faces to restore calculates a consistent set of
            //  saved faces (indices into savedFaces_) and saved vertices
            //  (indices into savedPoints_) to restore. The undoFaces have to
            //  be synced.
            void getUnrefimentSet
            (
                const labelList& undoFaces,
                labelList& localFaces,
                labelList& localPoints
            ) const;

            //- Restore selected faces and vertices.
            void setUnrefinement
            (
                const labelList& localFaces,
                const labelList& localPoints,
                polyTopoChange&
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
