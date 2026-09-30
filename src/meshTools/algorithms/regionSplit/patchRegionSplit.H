/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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
    Foam::patchRegionSplit

Description
    ...

SourceFiles
    patchRegionSplit.C

\*---------------------------------------------------------------------------*/

#ifndef patchRegionSplit_H
#define patchRegionSplit_H

#include "regionSplitBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;
class polyPatch;

/*---------------------------------------------------------------------------*\
                         Class patchRegionSplit Declaration
\*---------------------------------------------------------------------------*/

class patchRegionSplit
:
    public regionSplitBase
{
    // Private constructors

        //- Construct from a poly-mesh and a primitive patch
        template<class Patch>
        patchRegionSplit
        (
            const polyMesh& mesh,
            const Patch& patch
        );


public:

    //- Runtime type information
    ClassName("patchRegionSplit");


    // Constructors

        //- Construct from a poly-mesh and a set of faces
        patchRegionSplit
        (
            const polyMesh& mesh,
            const labelList& faces
        );

        //- Construct from a poly-patch
        patchRegionSplit(const polyPatch& patch);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
