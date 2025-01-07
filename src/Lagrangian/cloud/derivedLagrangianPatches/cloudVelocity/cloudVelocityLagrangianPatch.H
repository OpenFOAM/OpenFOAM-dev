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
    Foam::cloudVelocityLagrangianPatch

Description
    A cloud velocity Lagrangian patch. This patch defers to the cloud velocity
    boundary condition to determine the state of a particle after hitting a
    (non-constraint) patch. This allows the velocity boundary condition to
    control the rebound, escape, stick, etc..., behaviour. This is selected
    automatically by the base cloud. The user should never have to specify this
    patch explicitly.

SourceFiles
    cloudVelocityLagrangianPatch.C

\*---------------------------------------------------------------------------*/

#ifndef cloudVelocityLagrangianPatch_H
#define cloudVelocityLagrangianPatch_H

#include "LagrangianPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class cloudVelocityLagrangianPatch Declaration
\*---------------------------------------------------------------------------*/

class cloudVelocityLagrangianPatch
:
    public LagrangianPatch
{
public:

    //- Runtime type information
    TypeName("cloudVelocity");


    // Constructors

        //- Construct from a patch and a boundary mesh
        cloudVelocityLagrangianPatch
        (
            const polyPatch&,
            const LagrangianBoundaryMesh&
        );


    //- Destructor
    virtual ~cloudVelocityLagrangianPatch();


    // Member Functions

        //- Evaluate changes in elements that have tracked to this patch
        virtual void evaluate
        (
            PstreamBuffers&,
            LagrangianMesh&,
            const LagrangianScalarInternalDynamicField& fraction
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
