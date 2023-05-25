/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::filmSurfaceFvPatch

Description
    Foam::filmSurfaceFvPatch

SourceFiles
    filmSurfaceFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef filmSurfaceFvPatch_H
#define filmSurfaceFvPatch_H

#include "filmFvPatch.H"
#include "filmSurfacePolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class filmSurfaceFvPatch Declaration
\*---------------------------------------------------------------------------*/

class filmSurfaceFvPatch
:
    public filmFvPatch
{
    // Private member data

        mutable scalarField deltaCoeffs_;


public:

    //- Runtime type information
    TypeName(filmSurfacePolyPatch::typeName_());


    // Constructors

        //- Construct from components
        filmSurfaceFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm)
        :
            filmFvPatch(patch, bm)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
