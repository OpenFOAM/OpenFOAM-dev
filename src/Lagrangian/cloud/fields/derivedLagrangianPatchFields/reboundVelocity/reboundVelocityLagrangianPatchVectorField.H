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
    Foam::reboundVelocityLagrangianPatchVectorField

Description
    A rebound velocity boundary condition. Particles that hit a patch with this
    condition will rebound back into the volume of the mesh. The change in
    velocity is characterised by a coefficient of restitution and a fractional
    tangential momentum loss.

Usage
    \table
        Property | Description                         | Required? | Default
        e        | Coefficient of restitution          | yes       |
        mu       | Fractional tangential momentum loss | yes       |
    \endtable

    Example specification:
    \verbatim
    <LagrangianPatchName>
    {
        type            reboundVelocity;
        e               0.8;
        mu              10 [%];
    }
    \endverbatim

SourceFiles
    reboundVelocityLagrangianPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef reboundVelocityLagrangianPatchVectorField_H
#define reboundVelocityLagrangianPatchVectorField_H

#include "cloudVelocityLagrangianPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class reboundVelocityLagrangianPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class reboundVelocityLagrangianPatchVectorField
:
    public cloudVelocityLagrangianPatchVectorField
{
    // Private data

        //- Coefficient of restitution
        const scalar e_;

        //- Fractional tangential momentum loss
        const scalar mu_;


public:

    //- Runtime type information
    TypeName("reboundVelocity");


    // Constructors

        //- Construct from patch, internal field and dictionary
        reboundVelocityLagrangianPatchVectorField
        (
            const LagrangianPatch&,
            const regIOobject&,
            const dictionary&
        );

        //- Construct as copy
        reboundVelocityLagrangianPatchVectorField
        (
            const reboundVelocityLagrangianPatchVectorField&
        );

        //- Construct as copy setting the internal field reference
        reboundVelocityLagrangianPatchVectorField
        (
            const reboundVelocityLagrangianPatchVectorField&,
            const regIOobject&
        );

        //- Construct and return a clone
        virtual autoPtr<LagrangianPatchVectorField> clone() const
        {
            return autoPtr<LagrangianPatchVectorField>
            (
                new reboundVelocityLagrangianPatchVectorField(*this)
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
                new reboundVelocityLagrangianPatchVectorField
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

        //- Write
        virtual void write(Ostream&) const;

        //- Return the state after interaction with this patch
        virtual LagrangianState state() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
