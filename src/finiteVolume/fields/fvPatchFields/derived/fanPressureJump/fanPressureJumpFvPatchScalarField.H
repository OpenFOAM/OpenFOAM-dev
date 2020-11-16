/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::fanPressureJumpFvPatchScalarField

Description
    This boundary condition provides a pressure jump condition, using
    the \c cyclic condition as a base.

    The jump is specified as a \c Function1 type, to enable the use of, e.g.
    constant, polynomial, table values. This boundary condition can operate
    in two modes - standard and backward compatibility.
    In standard mode, jump is specified as a function of total volumetric
    flow rate through the patch. In backward compatibility mode, the boundary
    conditions serves as a direct replacement for the fanFvPatchField, where
    jump is defined as a function of local velocity.

Usage
    \table
        Property     | Description             | Required    | Default value
        patchType    | underlying patch type (should be \c cyclic) | yes |
        fanCurve     | fan curve function      | yes         |
        jumpTable    | jump table function (backward compatibility mode) | no |
        reverse      | reverse jump direction  | no          | false
        phi          | flux field name         | no          | phi
        rho          | density field name      | no          | rho
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            fanPressureJump;
        patchType       cyclic;
        fanCurve        table;
        fanCurveCoeffs
        {
            file            "$FOAM_CASE/constant/pressureVsQ";
            format          csv;
            nHeaderLine     1;
            refColumn       0;
            componentColumns 1(1);
            separator       ",";
            mergeSeparators no;
            outOfBounds     clamp;
            interpolationScheme linear;
        }
        value           uniform 0;
    }
    \endverbatim

    The above example shows the use of a comma separated (CSV) file to specify
    the jump condition.

Note
     The underlying \c patchType should be set to \c cyclic

See also
    Foam::Function1s

SourceFiles
    fanPressureJumpFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef fanPressureJumpFvPatchScalarField_H
#define fanPressureJumpFvPatchScalarField_H

#include "fixedJumpFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class fanPressureJumpFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class fanPressureJumpFvPatchScalarField
:
    public fixedJumpFvPatchScalarField
{
    // Private Data

        //- Fan curve
        autoPtr<Function1<scalar>> fanCurve_;

        //- Jump table (backward compatibility)
        autoPtr<Function1<scalar>> jumpTable_;

        //- Set true to reverse jump direction
        Switch reverse_;

        //- Name of the flux transporting the field
        word phiName_;

        //- Name of the density field used to normalise the mass flux
        //  if necessary
        word rhoName_;


    // Private Member Functions

        //- Calculate the fan pressure jump
        void calcFanJump();


public:

    //- Runtime type information
    TypeName("fanPressureJump");


    // Constructors

        //- Construct from patch and internal field
        fanPressureJumpFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fanPressureJumpFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patch onto a new patch
        fanPressureJumpFvPatchScalarField
        (
            const fanPressureJumpFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        fanPressureJumpFvPatchScalarField
        (
            const fanPressureJumpFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        fanPressureJumpFvPatchScalarField
        (
            const fanPressureJumpFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<scalar>>
            (
                new fanPressureJumpFvPatchScalarField(*this, iF)
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
