/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::mappedVelocityFluxFixedValueFvPatchField

Group
    grpInletBoundaryConditions grpCoupledBoundaryConditions

Description
    This boundary condition maps the velocity and flux from a neighbour patch
    to this patch

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        phi          | flux field name         | no          | phi
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            mappedVelocityFlux;
        phi             phi;
        value           uniform 0;  // place holder
    }
    \endverbatim

    The underlying sample mode should be set to \c nearestPatchFace or
    \c nearestFace

Note
    This boundary condition can only be applied to patches that are of
    the \c mappedPolyPatch type.

SeeAlso
    Foam::mappedPatchBase
    Foam::mappedPolyPatch
    Foam::mappedFvPatch
    Foam::fixedValueFvPatchVectorField

SourceFiles
    mappedVelocityFluxFixedValueFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef mappedVelocityFluxFixedValueFvPatchField_H
#define mappedVelocityFluxFixedValueFvPatchField_H

#include "fixedValueFvPatchFields.H"
#include "mappedFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class mappedVelocityFluxFixedValueFvPatch Declaration
\*---------------------------------------------------------------------------*/

class mappedVelocityFluxFixedValueFvPatchField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Name of flux field
        word phiName_;

public:

    //- Runtime type information
    TypeName("mappedVelocityFlux");


    // Constructors

        //- Construct from patch and internal field
        mappedVelocityFluxFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        mappedVelocityFluxFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  mappedVelocityFluxFixedValueFvPatchField
        //  onto a new patch
        mappedVelocityFluxFixedValueFvPatchField
        (
            const mappedVelocityFluxFixedValueFvPatchField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        mappedVelocityFluxFixedValueFvPatchField
        (
            const mappedVelocityFluxFixedValueFvPatchField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new mappedVelocityFluxFixedValueFvPatchField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        mappedVelocityFluxFixedValueFvPatchField
        (
            const mappedVelocityFluxFixedValueFvPatchField&,
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
                new mappedVelocityFluxFixedValueFvPatchField(*this, iF)
            );
        }


    // Member functions

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
