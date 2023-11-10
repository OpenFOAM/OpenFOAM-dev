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
    Foam::mappedVelocityFluxFvPatchField

Description
    This boundary condition maps the velocity and flux from a neighbouring
    patch to this patch.

Usage
    \table
        Property     | Description             | Required    | Default value
        phi          | flux field name         | no          | phi
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            mappedVelocityFlux;
        value           uniform (0 0 0);
    }
    \endverbatim

    This boundary condition can only be applied to patches that are of
    mappedPatchBase type.

See also
    Foam::mappedPatchBase
    Foam::mappedPolyPatch
    Foam::mappedFvPatch
    Foam::fixedValueFvPatchVectorField

SourceFiles
    mappedVelocityFluxFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef mappedVelocityFluxFvPatchField_H
#define mappedVelocityFluxFvPatchField_H

#include "fixedValueFvPatchFields.H"
#include "mappedFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class mappedVelocityFluxFvPatch Declaration
\*---------------------------------------------------------------------------*/

class mappedVelocityFluxFvPatchField
:
    public fixedValueFvPatchVectorField
{
    // Private Data

        //- Name of flux field
        word phiName_;


public:

    //- Runtime type information
    TypeName("mappedVelocityFlux");


    // Constructors

        //- Construct from patch, internal field and dictionary
        mappedVelocityFluxFvPatchField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  mappedVelocityFluxFvPatchField
        //  onto a new patch
        mappedVelocityFluxFvPatchField
        (
            const mappedVelocityFluxFvPatchField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        mappedVelocityFluxFvPatchField
        (
            const mappedVelocityFluxFvPatchField&
        ) = delete;

        //- Copy constructor setting internal field reference
        mappedVelocityFluxFvPatchField
        (
            const mappedVelocityFluxFvPatchField&,
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
                new mappedVelocityFluxFvPatchField(*this, iF)
            );
        }


    // Member Functions

        // Evaluation functions

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
