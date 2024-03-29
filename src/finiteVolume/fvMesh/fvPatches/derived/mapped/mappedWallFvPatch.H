/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::mappedWallFvPatch

Description
    Wall fv patch which can do interpolative mapping of values from another
    globally conforming fv patch

SourceFiles
    mappedWallFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef mappedWallFvPatch_H
#define mappedWallFvPatch_H

#include "wallFvPatch.H"
#include "mappedFvPatchBase.H"
#include "mappedWallPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class mappedWallFvPatch Declaration
\*---------------------------------------------------------------------------*/

class mappedWallFvPatch
:
    public wallFvPatch,
    public mappedFvPatchBase
{

public:

    //- Runtime type information
    TypeName(mappedWallPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        mappedWallFvPatch
        (
            const polyPatch& patch,
            const fvBoundaryMesh& bm
        )
        :
            wallFvPatch(patch, bm),
            mappedFvPatchBase(static_cast<const fvPatch&>(*this))
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
