/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::genericLagrangianPatch

Description
    Substitute for unknown patches. Used for post-processing when only
    basic LagrangianPatch info is needed.

SourceFiles
    genericLagrangianPatch.C

\*---------------------------------------------------------------------------*/

#ifndef genericLagrangianPatch_H
#define genericLagrangianPatch_H

#include "LagrangianPatch.H"
#include "genericPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class genericLagrangianPatch Declaration
\*---------------------------------------------------------------------------*/

class genericLagrangianPatch
:
    public LagrangianPatch
{
public:

    //- Runtime type information
    TypeName(genericPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        genericLagrangianPatch
        (
            const polyPatch& patch,
            const LagrangianBoundaryMesh& bm
        )
        :
            LagrangianPatch(patch, bm)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
