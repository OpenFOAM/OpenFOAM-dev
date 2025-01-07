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
    Foam::stickVelocityLagrangianPatchVectorField

Description
    A stick velocity boundary condition. Particles that hit a patch with this
    condition will stick to it and remain at the patch's velocity (stationary
    if the mesh is not moving) from then on.

Usage
    Example specification:
    \verbatim
    <LagrangianPatchName>
    {
        type            stickVelocity;
    }
    \endverbatim

SourceFiles
    stickVelocityLagrangianPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef stickVelocityLagrangianPatchVectorField_H
#define stickVelocityLagrangianPatchVectorField_H

#include "cloudVelocityLagrangianPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class stickVelocityLagrangianPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class stickVelocityLagrangianPatchVectorField
:
    public cloudVelocityLagrangianPatchVectorField
{
public:

    //- Runtime type information
    TypeName("stickVelocity");


    // Constructors

        //- Inherit constructors
        using
            cloudVelocityLagrangianPatchVectorField::
            cloudVelocityLagrangianPatchVectorField;

        //- Construct and return a clone
        virtual autoPtr<LagrangianPatchVectorField> clone() const
        {
            return autoPtr<LagrangianPatchVectorField>
            (
                new stickVelocityLagrangianPatchVectorField(*this)
            );
        }

        //- Construct and return a clone setting the internal field reference
        virtual autoPtr<LagrangianPatchVectorField> clone
        (
            const regIOobject& iIo
        ) const
        {
            return autoPtr<LagrangianPatchVectorField>
            (
                new stickVelocityLagrangianPatchVectorField
                (
                    *this,
                    iIo
                )
            );
        }


    // Member Functions

        //- Evaluate the patch field
        virtual void evaluate
        (
            PstreamBuffers&,
            const LagrangianScalarInternalDynamicField& fraction
        );

        //- Return the state after interaction with this patch
        virtual LagrangianState state() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
