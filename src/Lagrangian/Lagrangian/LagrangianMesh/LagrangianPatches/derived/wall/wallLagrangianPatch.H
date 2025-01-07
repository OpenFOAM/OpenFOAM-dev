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
    Foam::wallLagrangianPatch

Description
    Wall Lagrangian patch. Provides no additional functionality over the base
    patch. Is simply a type that can be identified as a wall.

SourceFiles
    wallLagrangianPatch.C

\*---------------------------------------------------------------------------*/

#ifndef wallLagrangianPatch_H
#define wallLagrangianPatch_H

#include "LagrangianPatch.H"
#include "wallPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class wallLagrangianPatch Declaration
\*---------------------------------------------------------------------------*/

class wallLagrangianPatch
:
    public LagrangianPatch
{
public:

    //- Runtime type information
    TypeName(wallPolyPatch::typeName_());


    // Constructors

        //- Construct from a patch and a boundary mesh
        wallLagrangianPatch(const polyPatch&, const LagrangianBoundaryMesh&);


    //- Destructor
    virtual ~wallLagrangianPatch();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
