/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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
    Foam::mappedExtrudedWallPointPatch

SourceFiles
    mappedExtrudedWallPointPatch.C

\*---------------------------------------------------------------------------*/

#ifndef mappedExtrudedWallPointPatch_H
#define mappedExtrudedWallPointPatch_H

#include "wallPointPatch.H"
#include "mappedExtrudedWallPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class mappedExtrudedWallPointPatch Declaration
\*---------------------------------------------------------------------------*/

class mappedExtrudedWallPointPatch
:
    public wallPointPatch
{

public:

    //- Runtime type information
    TypeName(mappedExtrudedWallPolyPatch::typeName_());


    // Constructors

        //- Construct from polyPatch
        mappedExtrudedWallPointPatch
        (
            const polyPatch& patch,
            const pointBoundaryMesh& bm
        )
        :
            wallPointPatch(patch, bm)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
