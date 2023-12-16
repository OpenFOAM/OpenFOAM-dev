/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::fvBoundaryMesh

Description
    Foam::fvBoundaryMesh

SourceFiles
    fvBoundaryMesh.C

\*---------------------------------------------------------------------------*/

#ifndef fvBoundaryMesh_H
#define fvBoundaryMesh_H

#include "fvPatchList.H"
#include "lduInterfacePtrsList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;
class polyBoundaryMesh;

/*---------------------------------------------------------------------------*\
                       Class fvBoundaryMesh Declaration
\*---------------------------------------------------------------------------*/

class fvBoundaryMesh
:
    public fvPatchList
{
    // Private Data

        //- Reference to mesh
        const fvMesh& mesh_;


    // Private Member Functions

        //- Add fvPatches corresponding to the given polyBoundaryMesh
        void addPatches(const polyBoundaryMesh&);


protected:

        //- Update boundary based on new polyBoundaryMesh
        void readUpdate(const polyBoundaryMesh&);


public:

    friend class fvMesh;


    // Constructors

        //- Construct with zero size
        fvBoundaryMesh(const fvMesh&);

        //- Construct from polyBoundaryMesh
        fvBoundaryMesh(const fvMesh&, const polyBoundaryMesh&);

        //- Move constructor
        fvBoundaryMesh(fvBoundaryMesh&&) = default;

        //- Disallow default bitwise copy construction
        fvBoundaryMesh(const fvBoundaryMesh&) = delete;


    // Member Functions

        // Access

            //- Return the mesh reference
            const fvMesh& mesh() const
            {
                return mesh_;
            }

            //- Return a list of pointers for each patch
            //  with only those pointing to interfaces being set
            lduInterfacePtrsList interfaces() const;

            //- Find patch index given a name
            label findIndex(const word& patchName) const;

            //- Find patch indices given a name
            labelList findIndices(const wordRe&, const bool useGroups) const;


        //- Correct patches after moving points
        void movePoints();

        //- Reorders patches. Ordering does not have to be done in
        //  ascending or descending order. Reordering has to be unique.
        //  (is shuffle). If validBoundary call needs to be parallel
        //  sync
        void shuffle(const labelUList& newToOld, const bool validBoundary);


    // Member Operators

        //- Return const and non-const reference to fvPatch by index.
        using fvPatchList::operator[];

        //- Return const reference to fvPatch by name.
        const fvPatch& operator[](const word&) const;

        //- Return reference to fvPatch by name.
        fvPatch& operator[](const word&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fvBoundaryMesh&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
