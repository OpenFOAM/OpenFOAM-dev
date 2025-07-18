/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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
        Property     | Description             | Required    | Default value
        phi          | flux field name         | no          | phi
        tangentialVelocity | tangential velocity | no  |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            pressureInletOutletVelocity;
        phi             phi;
        tangentialVelocity (0 0 0);
        value           uniform (0 0 0);
    }
    \endverbatim

    The \c tangentialVelocity entry is a Function1 of time, see
    Foam::Function1s.

See also
    Foam::directionMixedFvPatchVectorField
    Foam::totalPressureFvPatchScalarField
    Foam::entrainmentPressureFvPatchScalarField
    Foam::fixedValueFvPatchField

SourceFiles
    pressureInletOutletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef pressureInletOutletVelocityFvPatchVectorField_H
#define pressureInletOutletVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "directionMixedFvPatchFields.H"
#include "Function1.H"

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

        //- Flux field name
        word phiName_;

        //- Optional tangential velocity component
        autoPtr<Function1<vector>> tangentialVelocity_;


public:

    //- Runtime type information
    TypeName("pressureInletOutletVelocity");


    // Constructors

        //- Construct from patch, internal field and dictionary
        pressureInletOutletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  pressureInletOutletVelocityFvPatchVectorField onto a new patch
        pressureInletOutletVelocityFvPatchVectorField
        (
            const pressureInletOutletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
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
            const autoPtr<Function1<vector>>& tangentialVelocity() const
            {
                return tangentialVelocity_;
            }


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member Operators

        virtual void operator=(const fvPatchField<vector>& pvf);

        //- Inherit assignment
        using directionMixedFvPatchVectorField::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
