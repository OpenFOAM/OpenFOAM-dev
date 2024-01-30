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
    Foam::mappedFvPatch

Description
    Fv patch which can do interpolative mapping of values from another globally
    conforming fv patch

SourceFiles
    mappedFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef mappedFvPatch_H
#define mappedFvPatch_H

#include "fvPatch.H"
#include "mappedFvPatchBase.H"
#include "mappedPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class mappedFvPatch Declaration
\*---------------------------------------------------------------------------*/

class mappedFvPatch
:
    public fvPatch,
    public mappedFvPatchBase
{

public:

    //- Runtime type information
    TypeName(mappedPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        mappedFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm)
        :
            fvPatch(patch, bm),
            mappedFvPatchBase(static_cast<const fvPatch&>(*this))
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
