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
    Foam::nonConformalMappedWallFvPatch

Description
    Wall fv patch which can do non-conformal mapping of values from another
    potentially non-globally conforming wall fv patch. As nonConformalFvPatch,
    but the neighbouring patch is local and known and is made available by this
    class.

See also
    Foam::nonConformalFvPatch
    Foam::nonConformalMappedPolyPatch

SourceFiles
    nonConformalMappedWallFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalMappedWallFvPatch_H
#define nonConformalMappedWallFvPatch_H

#include "wallFvPatch.H"
#include "nonConformalFvPatch.H"
#include "nonConformalMappedFvPatchBase.H"
#include "nonConformalMappedWallPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class nonConformalMappedWallFvPatch Declaration
\*---------------------------------------------------------------------------*/

class nonConformalMappedWallFvPatch
:
    public wallFvPatch,
    public nonConformalFvPatch,
    public nonConformalMappedFvPatchBase
{
    // Private Data

        //- Reference to the polyPatch
        const nonConformalMappedWallPolyPatch& nonConformalMappedWallPolyPatch_;


public:

    //- Runtime type information
    TypeName(nonConformalMappedWallPolyPatch::typeName_());


    // Constructors

        //- Construct from polyPatch and fvBoundaryMesh
        nonConformalMappedWallFvPatch
        (
            const polyPatch& patch,
            const fvBoundaryMesh& bm
        );


    //- Destructor
    virtual ~nonConformalMappedWallFvPatch();


    // Member Functions

        // Access

            //- Inherit wallFvPatch patch method
            using wallFvPatch::patch;

            //- Poly patch
            const nonConformalMappedWallPolyPatch&
                nonConformalMappedWallPatch() const;

            //- Neighbour patch
            const nonConformalMappedWallFvPatch& nbrPatch() const;

            //- Is this patch the owner?
            bool owner() const;

            //- Return face face-poly-faces
            const labelList& polyFaces() const;

            //- Return the start label of this patch in the polyMesh face list.
            //  Raises an error unless the patch is empty.
            virtual label start() const;

            //- Return the size
            virtual label size() const;

            //- Return the face-cells
            virtual const labelUList& faceCells() const;

            //- Return the face-poly-faces patch field type
            virtual word polyFacesType() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
