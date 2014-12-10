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
    Foam::meshTriangulation

Description
    Triangulation of mesh faces. Generates (multiply connected) trisurface.

    All patch faces keep their patchID as triangle region.
    Internal faces get the supplied region number.

SourceFiles
    meshTriangulation.C

\*---------------------------------------------------------------------------*/


#ifndef meshTriangulation_H
#define meshTriangulation_H

#include "triSurface.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class primitiveMesh;

/*---------------------------------------------------------------------------*\
                           Class meshTriangulation Declaration
\*---------------------------------------------------------------------------*/

class meshTriangulation
:
    public triSurface
{
    // Private data

        //- Number of triangles in this that are internal to the surface.
        label nInternalFaces_;

        //- From triangle to mesh face
        labelList faceMap_;

    // Private Member Functions

        //- Is face internal to the subset.
        static bool isInternalFace
        (
            const primitiveMesh&,
            const boolList& includedCell,
            const label faceI
        );

        //- find boundary faces of subset.
        static void getFaces
        (
            const primitiveMesh&,
            const boolList& includedCell,
            boolList& faceIsCut,
            label& nFaces,
            label& nInternalFaces
        );

        //- Add triangulation of face to triangles. Optionally reverse.
        void insertTriangles
        (
            const triFaceList&,
            const label faceI,
            const label regionI,
            const bool reverse,

            List<labelledTri>& triangles,
            label& triI
        );


public:

    ClassName("meshTriangulation");


    // Constructors

        //- Construct null
        meshTriangulation();

        //- Construct from selected mesh cell and region number to be used
        //  for triangles resulting from internal faces. (all boundary triangles
        //  get polyMesh patch id).
        //  faceCentreDecomposition = true  : decomposition around face centre
        //                            false : decomposition using
        //                                    existing vertices
        meshTriangulation
        (
            const polyMesh&,
            const label internalFacesPatch,
            const boolList& includedCell,
            const bool faceCentreDecomposition = false
        );


    // Member Functions

        //- number of triangles in *this which are internal to the surface
        label nInternalFaces() const
        {
            return nInternalFaces_;
        }

        //- from triangle to mesh face
        const labelList& faceMap() const
        {
            return faceMap_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
