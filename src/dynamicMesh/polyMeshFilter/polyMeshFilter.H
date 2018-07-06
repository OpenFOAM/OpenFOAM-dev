/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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
    Foam::polyMeshFilter

Description
    Remove the edges and faces of a polyMesh whilst satisfying the given mesh
    quality criteria.

    Works on a copy of the mesh.

SourceFiles
    polyMeshFilter.C
    polyMeshFilterTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef polyMeshFilter_H
#define polyMeshFilter_H

#include "IOdictionary.H"
#include "Time.H"
#include "List.H"
#include "autoPtr.H"
#include "scalarField.H"
#include "polyMeshFilterSettings.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;
class fvMesh;
class PackedBoolList;
class faceSet;

/*---------------------------------------------------------------------------*\
                       Class polyMeshFilter Declaration
\*---------------------------------------------------------------------------*/

class polyMeshFilter
:
    private polyMeshFilterSettings
{
    // Private data

        //- Reference to the original mesh
        const fvMesh& mesh_;

        //- Copy of the original mesh to perform the filtering on
        autoPtr<fvMesh> newMeshPtr_;

        //- Original point priorities. If a point has a higher priority than
        //  another point then the edge between them collapses towards the
        //  point with the higher priority. (e.g. 2 is a higher priority than 1)
        labelList originalPointPriority_;

        //- Point priority associated with the new mesh
        autoPtr<labelList> pointPriority_;

        //- The minimum edge length for each edge
        scalarField minEdgeLen_;

        //- The face filter factor for each face
        scalarField faceFilterFactor_;


    // Private Member Functions

        template<class T>
        static void updateSets(const mapPolyMesh& map);

        template<class T>
        static void copySets(const polyMesh& oldMesh, const polyMesh& newMesh);

        label filterFacesLoop(const label nOriginalBadFaces);

        label filterFaces
        (
            polyMesh& newMesh,
            scalarField& newMeshFaceFilterFactor,
            labelList& origToCurrentPointMap
        );

        label filterEdges
        (
            polyMesh& newMesh,
            scalarField& newMeshMinEdgeLen,
            labelList& origToCurrentPointMap
        );

        //- Increment pointErrorCount for points attached to a bad face
        void updatePointErrorCount
        (
            const PackedBoolList& isErrorPoint,
            const labelList& oldToNewMesh,
            labelList& pointErrorCount
        ) const;


        //- Given the new points that are part of bad faces, and a map from the
        //  old mesh points to the new mesh points, relax minEdgeLen_
        void checkMeshEdgesAndRelaxEdges
        (
            const polyMesh& newMesh,
            const labelList& oldToNewMesh,
            const PackedBoolList& isErrorPoint,
            const labelList& pointErrorCount
        );

        //- Given the new points that are part of bad faces, and a map from the
        //  old mesh points to the new mesh points, relax faceFilterFactor_
        void checkMeshFacesAndRelaxEdges
        (
            const polyMesh& newMesh,
            const labelList& oldToNewMesh,
            const PackedBoolList& isErrorPoint,
            const labelList& pointErrorCount
        );

        // Mark boundary points
        // boundaryPoint:
        // + -1 : point not on boundary
        // +  0 : point on a real boundary
        // + >0 : point on a processor patch with that ID
        // TODO: Need to mark boundaryEdges as well, as an edge may have two
        //       boundary points but not itself lie on a boundary
        void updatePointPriorities
        (
            const polyMesh& newMesh,
            const labelList& pointMap
        );

        //- Print min/mean/max data for a field
        void printScalarFieldStats
        (
            const string desc,
            const scalarField& fld
        ) const;

        //- Update minEdgeLen_ for the new mesh based upon the movement of the
        //  old points to the new points
        void mapOldMeshEdgeFieldToNewMesh
        (
            const polyMesh& newMesh,
            const labelList& pointMap,
            scalarField& newMeshMinEdgeLen
        ) const;

        //- Update faceFilterFactor_ for the new mesh based upon the movement
        //  of the old faces to the new faces
        void mapOldMeshFaceFieldToNewMesh
        (
            const polyMesh& newMesh,
            const labelList& faceMap,
            scalarField& newMeshFaceFilterFactor
        ) const;

        //- Maintain a map of the original mesh points to the latest version of
        //  the filtered mesh.
        void updateOldToNewPointMap
        (
            const labelList& currToNew,
            labelList& origToCurrentPointMap
        ) const;

        //- Disallow default bitwise copy construct
        polyMeshFilter(const polyMeshFilter&);

        //- Disallow default bitwise assignment
        void operator=(const polyMeshFilter&);


public:

    //- Runtime type information
    ClassName("polyMeshFilter");


    // Constructors

        //- Construct from fvMesh
        explicit polyMeshFilter(const fvMesh& mesh);

        //- Construct from fvMesh and a label list of point priorities
        polyMeshFilter(const fvMesh& mesh, const labelList& pointPriority);

        //- Construct from fvMesh and a label list of point priorities
        polyMeshFilter
        (
            const fvMesh& mesh,
            const labelList& pointPriority,
            const dictionary& dict
        );


    //- Destructor
    ~polyMeshFilter();


    // Member Functions

        // Access

            //- Return reference to the filtered mesh. Does not check if the
            //  mesh has actually been filtered.
            const autoPtr<fvMesh>& filteredMesh() const;

            //- Return the new pointPriority list.
            const autoPtr<labelList>& pointPriority() const;


        // Edit

            //- Return a copy of an fvMesh
            static autoPtr<fvMesh> copyMesh(const fvMesh& mesh);

            //- Copy loaded topoSets from the old mesh to the new mesh
            static void copySets
            (
                const polyMesh& oldMesh,
                const polyMesh& newMesh
            );

            //- Update the loaded topoSets
            static void updateSets(const mapPolyMesh& map);

            //- Filter edges and faces
            label filter(const label nOriginalBadFaces);

            //- Filter all faces that are in the face set
            label filter(const faceSet& fSet);

            //- Filter edges only.
            label filterEdges(const label nOriginalBadFaces);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "polyMeshFilterTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
