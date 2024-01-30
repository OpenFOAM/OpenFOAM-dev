/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
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
    Foam::nonConformalErrorFvPatch

Description
    Non-conformal error FV patch. As nonConformalFvPatch. This patch is a
    non-coupled non-conformal patch which is used to manage the errors created
    during the construction of a non-conformal coupled interface. Every patch
    used as the original patch of the owner side of a non-conformal coupled
    interface must also have an associated error patch.

See also
    Foam::nonConformalFvPatch

SourceFiles
    nonConformalErrorFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalErrorFvPatch_H
#define nonConformalErrorFvPatch_H

#include "fvPatch.H"
#include "nonConformalErrorPolyPatch.H"
#include "nonConformalFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class nonConformalErrorFvPatch Declaration
\*---------------------------------------------------------------------------*/

class nonConformalErrorFvPatch
:
    public fvPatch,
    public nonConformalFvPatch
{
    // Private Data

        //- Reference to the polyPatch
        const nonConformalErrorPolyPatch& nonConformalErrorPolyPatch_;


public:

    //- Runtime type information
    TypeName(nonConformalErrorPolyPatch::typeName_());


    // Constructors

        //- Construct from polyPatch and fvBoundaryMesh
        nonConformalErrorFvPatch
        (
            const polyPatch& patch,
            const fvBoundaryMesh& bm
        );


    //- Destructor
    virtual ~nonConformalErrorFvPatch();


    // Member Functions

        // Access

            //- Inherit fvPatch patch method
            using fvPatch::patch;

            //- Poly patch
            const nonConformalErrorPolyPatch& nonConformalErrorPatch() const;

            //- Return face face-poly-faces
            const labelList& polyFaces() const;

            //- Return the start label of this patch in the polyMesh face list.
            //  Raises an error unless the patch is empty.
            virtual label start() const;

            //- Return the size
            virtual label size() const;

            //- Return the face-cells
            virtual const labelUList& faceCells() const;

            //- Return delta vectors. It doesn't really matter what these are,
            //  this is an error patch after all, they just have to be
            //  stable/non-zero. The standard formulation in fvPatch is not
            //  stable because it depends on the face normal, which for an
            //  error face is arbitrary. So, we use the same as for coupled
            //  patches; the vector between the cell and the face centres.
            virtual tmp<vectorField> delta() const;

            //- Return the face-poly-faces patch field type
            virtual word polyFacesType() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
