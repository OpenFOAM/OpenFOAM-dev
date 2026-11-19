/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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
    Foam::pressureInletOutletVelocityFvPatchVectorField

Description
    Velocity inlet/outlet boundary condition for patches where the pressure is
    specified in some manner, e.g. fixedValue, totalPressure,
    entrainmentPressure etc.

    This boundary condition is directionMixed:
        Flux outflow: zeroGradient condition
        Flux inflow: zeroGradient condition for normal component
                     fixedValue for tangential component

    The external tangential velocity can be optionally specified otherwise it
    is assumed zero.

Usage
    \table
        Property     | Description               | Required | Default value
        phi          | flux field name           | no       | phi
        tangentialVelocity | tangential velocity | no       |
        value        | initial velocity          | no       | evaluated
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            pressureInletOutletVelocity;
        tangentialVelocity $internalField;
    }
    \endverbatim

    The \c tangentialVelocity entry is a \c UFunctionalDimensionedField
    providing total flexibility in the specification of the tangential velocity,
    the normal component is removed automatically from the field returned by the
    function.

See also
    Foam::DimensionedFieldFunction
    Foam::UFunctionalDimensionedField
    Foam::directionMixedFvPatchVectorField
    Foam::totalPressureFvPatchScalarField
    Foam::entrainmentPressureFvPatchScalarField

SourceFiles
    pressureInletOutletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef pressureInletOutletVelocityFvPatchVectorField_H
#define pressureInletOutletVelocityFvPatchVectorField_H

#include "directionMixedFvPatchFields.H"
#include "UFunctionalDimensionedField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class pressureInletOutletVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class pressureInletOutletVelocityFvPatchVectorField
:
    public directionMixedFvPatchVectorField
{
    // Private Data

        //- Optional tangential velocity
        UFunctionalDimensionedField<vector, fvPatch> tangentialVelocity_;

        //- Flux field name
        word phiName_;


public:

    //- Runtime type information
    TypeName("pressureInletOutletVelocity");


    // Constructors

        //- Construct from patch, internal field and dictionary
        pressureInletOutletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, fvMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  pressureInletOutletVelocityFvPatchVectorField onto a new patch
        pressureInletOutletVelocityFvPatchVectorField
        (
            const pressureInletOutletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, fvMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        pressureInletOutletVelocityFvPatchVectorField
        (
            const pressureInletOutletVelocityFvPatchVectorField&
        ) = delete;

        //- Copy constructor setting internal field reference
        pressureInletOutletVelocityFvPatchVectorField
        (
            const pressureInletOutletVelocityFvPatchVectorField&,
            const DimensionedField<vector, fvMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, fvMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new pressureInletOutletVelocityFvPatchVectorField(*this, iF)
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

            //- Return the tangential velocity Function1
            const DimensionedField<vector, fvPatch>& tangentialVelocity() const
            {
                return tangentialVelocity_;
            }


        // Mapping functions

            //- Map the given fvPatchField onto this fvPatchField
            virtual void map(const fvPatchField<vector>&, const fieldMapper&);

            //- Reset the fvPatchField to the given fvPatchField
            //  Used for mesh to mesh mapping
            virtual void reset(const fvPatchField<vector>&);


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member Operators

        //- Virtual assign to an expression
        virtual void operator=(const expression::Assignment<vectorField>&);

        //- Inherit assignment
        using directionMixedFvPatchVectorField::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
