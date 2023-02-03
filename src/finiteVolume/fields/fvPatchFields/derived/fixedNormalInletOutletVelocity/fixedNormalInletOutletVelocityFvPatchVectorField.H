/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2023 OpenFOAM Foundation
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
    Foam::fixedNormalInletOutletVelocityFvPatchVectorField

Description

    This velocity inlet/outlet boundary condition combines a fixed normal
    component obtained from the "normalVelocity" patchField supplied with a
    fixed or zero-gradiented tangential component depending on the direction
    of the flow and the setting of "fixTangentialInflow":
    - Outflow: apply zero-gradient condition to tangential components
    - Inflow:
      - fixTangentialInflow is true
        apply value provided by the normalVelocity patchField to the
        tangential components
      - fixTangentialInflow is false
        apply zero-gradient condition to tangential components.

Usage
    \table
        Property     | Description             | Required    | Default value
        phi          | flux field name         | no          | phi
  fixTangentialInflow | If true fix the tangential component for inflow | yes |
        normalVelocity | patchField providing the normal velocity field | yes |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            fixedNormalInletOutletVelocity;

        fixTangentialInflow false;
        normalVelocity
        {
            type            uniformFixedValue;
            uniformValue    sine;
            uniformValueCoeffs
            {
                frequency 1;
                amplitude table
                (
                    (0  0)
                    (2  0.088)
                    (8  0.088)
                );
                scale     (0 1 0);
                level     (0 0 0);
            }
        }

        value           uniform (0 0 0);
    }
    \endverbatim

SourceFiles
    fixedNormalInletOutletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedNormalInletOutletVelocityFvPatchVectorField_H
#define fixedNormalInletOutletVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "directionMixedFvPatchFields.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class fixedNormalInletOutletVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class fixedNormalInletOutletVelocityFvPatchVectorField
:
    public directionMixedFvPatchVectorField
{
    // Private Data

        //- Flux field name
        word phiName_;

        //- Set true to fix the tangential component for inflow
        Switch fixTangentialInflow_;

        //- BC which provided the normal component of the velocity
        tmp<fvPatchVectorField> normalVelocity_;


public:

    //- Runtime type information
    TypeName("fixedNormalInletOutletVelocity");


    // Constructors

        //- Construct from patch and internal field
        fixedNormalInletOutletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedNormalInletOutletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  fixedNormalInletOutletVelocityFvPatchVectorField onto a new patch
        fixedNormalInletOutletVelocityFvPatchVectorField
        (
            const fixedNormalInletOutletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        fixedNormalInletOutletVelocityFvPatchVectorField
        (
            const fixedNormalInletOutletVelocityFvPatchVectorField&
        ) = delete;

        //- Copy constructor setting internal field reference
        fixedNormalInletOutletVelocityFvPatchVectorField
        (
            const fixedNormalInletOutletVelocityFvPatchVectorField&,
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
                new fixedNormalInletOutletVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member Functions

        // Attributes

            //- Return true: this patch field is altered by assignment
            virtual bool assignable() const
            {
                return true;
            }


        // Access

            //- Return the name of phi
            const word& phiName() const
            {
                return phiName_;
            }

            //- Return reference to the name of phi to allow adjustment
            word& phiName()
            {
                return phiName_;
            }

            Switch fixTangentialInflow() const
            {
                return fixTangentialInflow_;
            }

            //- Return the BC which provides the normal component of velocity
            const fvPatchVectorField& normalVelocity() const
            {
                return normalVelocity_();
            }


        // Mapping functions

            //- Map the given fvPatchField onto this fvPatchField
            virtual void map
            (
                const fvPatchVectorField&,
                const fvPatchFieldMapper&
            );

            //- Reset the fvPatchField to the given fvPatchField
            //  Used for mesh to mesh mapping
            virtual void reset(const fvPatchVectorField&);


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member Operators

        virtual void operator=(const fvPatchField<vector>& pvf);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
