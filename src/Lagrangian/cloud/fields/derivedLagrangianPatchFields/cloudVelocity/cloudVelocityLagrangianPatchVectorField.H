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
    Foam::cloudVelocityLagrangianPatchVectorField

Description
    Base class for cloud velocity conditions. Accessed by the cloud velocity
    Lagrangian patch. Provides an interface to return the state of a particle
    after hitting a (non-constraint) patch. This allows derivations to control
    the rebound, escape, stick, etc..., behaviour.

See also
    Foam::cloudVelocityLagrangianPatch

SourceFiles
    cloudVelocityLagrangianPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef cloudVelocityLagrangianPatchVectorField_H
#define cloudVelocityLagrangianPatchVectorField_H

#include "LagrangianPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class cloudVelocityLagrangianPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class cloudVelocityLagrangianPatchVectorField
:
    public LagrangianPatchVectorField
{
public:

    //- Runtime type information
    TypeName("cloudVelocity");


    // Constructors

        //- Inherit constructors
        using LagrangianPatchVectorField::LagrangianPatchVectorField;


    // Member Functions

        //- Return the state after interaction with this patch
        virtual LagrangianState state() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
