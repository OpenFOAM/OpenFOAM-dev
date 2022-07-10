/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::pointBoundaryMesh

Description
    Foam::pointBoundaryMesh

SourceFiles
    pointBoundaryMesh.C

\*---------------------------------------------------------------------------*/

#ifndef pointBoundaryMesh_H
#define pointBoundaryMesh_H

#include "pointPatchList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class pointMesh;
class polyBoundaryMesh;

/*---------------------------------------------------------------------------*\
                      Class pointBoundaryMesh Declaration
\*---------------------------------------------------------------------------*/

class pointBoundaryMesh
:
    public pointPatchList
{
    // Private Data

        //- Reference to mesh
        const pointMesh& mesh_;


    // Private Member Functions

        //- Calculate the geometry for the patches (transformation tensors etc.)
        void calcGeometry();


public:

    //- Declare friendship with pointMesh
    friend class pointMesh;


    // Constructors

        //- Construct from polyBoundaryMesh
        pointBoundaryMesh
        (
            const pointMesh&,
            const polyBoundaryMesh&
        );

        //- Disallow default bitwise copy construction
        pointBoundaryMesh(const pointBoundaryMesh&) = delete;


    // Member Functions

        //- Return the mesh reference
        const pointMesh& mesh() const
        {
            return mesh_;
        }

        //- Find patch index given a name
        label findPatchID(const word& patchName) const;

        //- Find patch indices given a name
        labelList findIndices(const wordRe&, const bool useGroups) const;

        //- Correct pointBoundaryMesh after moving points
        void movePoints(const pointField&);

        //- Correct pointBoundaryMesh after topology update
        void topoChange();

        //- Reset pointBoundaryMesh with respect to the updated polyBoundaryMesh
        //  For run-time mesh replacement and mesh to mesh mapping
        void reset();

        //- Reorders patches. Ordering does not have to be done in
        //  ascending or descending order. Reordering has to be unique.
        //  (is shuffle). If validBoundary calls topoChange()
        //  after reordering to recalculate data (so call needs to be parallel
        //  sync in that case)
        void shuffle(const labelUList& newToOld, const bool validBoundary);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const pointBoundaryMesh&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
