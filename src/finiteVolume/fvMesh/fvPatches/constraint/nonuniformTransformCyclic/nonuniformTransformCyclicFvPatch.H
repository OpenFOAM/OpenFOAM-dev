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
    Foam::nonuniformTransformCyclicFvPatch

Description
    Cyclic-plane patch.

SourceFiles
    nonuniformTransformCyclicFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef nonuniformTransformCyclicFvPatch_H
#define nonuniformTransformCyclicFvPatch_H

#include "cyclicFvPatch.H"
#include "nonuniformTransformCyclicPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class nonuniformTransformCyclicFvPatch Declaration
\*---------------------------------------------------------------------------*/

class nonuniformTransformCyclicFvPatch
:
    public cyclicFvPatch
{

public:

    //- Runtime type information
    TypeName(nonuniformTransformCyclicPolyPatch::typeName_());


    // Constructors

        //- Construct from polyPatch
        nonuniformTransformCyclicFvPatch
        (
            const polyPatch& patch,
            const fvBoundaryMesh& bm
        )
        :
            cyclicFvPatch(patch, bm)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
