/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::mappedFlowRateVelocityFvPatchVectorField

Description
    This boundary condition maps the flow rate from a neighbouring patch to
    this patch, then uses it to set a corresponding velocity in a direction
    normal to the patch.

Usage
    \table
        Property     | Description             | Required    | Default value
        nbrPhi       | name of flux field on neighbour mesh | no | phi
        phi          | flux field name         | no          | phi
        rho          | density field name      | no          | rho
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            mappedFlowRateVelocity;
        value           uniform (0 0 0);
    }
    \endverbatim

    This boundary condition can only be applied to patches that are of
    mapped type.

See also
    Foam::mappedPolyPatch
    Foam::mappedFvPatch
    Foam::fixedValueFvPatchVectorField

SourceFiles
    mappedFlowRateVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef mappedFlowRateVelocityFvPatchVectorField_H
#define mappedFlowRateVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
            Class mappedFlowRateVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class mappedFlowRateVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private Data

        //- Name of the neighbor flux setting the inlet mass flux
        word nbrPhiName_;

        //- Name of the local mass flux
        word phiName_;

        //- Name of the density field used to normalise the mass flux
        word rhoName_;


public:

   //- Runtime type information
   TypeName("mappedFlowRateVelocity");


   // Constructors

        //- Construct from patch, internal field and dictionary
        mappedFlowRateVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  mappedFlowRateVelocityFvPatchVectorField
        //  onto a new patch
        mappedFlowRateVelocityFvPatchVectorField
        (
            const mappedFlowRateVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        mappedFlowRateVelocityFvPatchVectorField
        (
            const mappedFlowRateVelocityFvPatchVectorField&
        ) = delete;

        //- Copy constructor setting internal field reference
        mappedFlowRateVelocityFvPatchVectorField
        (
            const mappedFlowRateVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new mappedFlowRateVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member Functions


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
