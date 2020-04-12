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
    Foam::compressible::alphatWallFunctionFvPatchScalarField

Description
    This boundary condition provides a turbulent thermal diffusivity condition
    when using wall functions
    - replicates OpenFOAM v1.5 (and earlier) behaviour

    The turbulent thermal diffusivity calculated using:

        \f[
            \alpha_t = \frac{\mu_t}{Pr_t}
        \f]

    where

    \vartable
        \alpha_t| turbulence thermal diffusivity
        \mu_t   | turbulence viscosity
        Pr_t    | turblent Prandtl number
    \endvartable

Usage
    \table
        Property     | Description             | Required    | Default value
        nut          | turbulence viscosity field name | no  | nut
        Prt          | turbulent Prandtl number | no          | 0.85
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            alphatWallFunction;
        nut             nut;
        Prt             0.85;
        value           uniform 0; // optional value entry
    }
    \endverbatim

See also
    Foam::fixedValueFvPatchField

SourceFiles
    alphatWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef alphatWallFunctionFvPatchScalarField_H
#define alphatWallFunctionFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
            Class alphatWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class alphatWallFunctionFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private Data

        //- Turbulent Prandtl number (default = 0.85)
        scalar Prt_;


public:

    //- Runtime type information
    TypeName("compressible::alphatWallFunction");


    // Constructors

        //- Construct from patch and internal field
        alphatWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        alphatWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  alphatWallFunctionFvPatchScalarField
        //  onto a new patch
        alphatWallFunctionFvPatchScalarField
        (
            const alphatWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        alphatWallFunctionFvPatchScalarField
        (
            const alphatWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new alphatWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        alphatWallFunctionFvPatchScalarField
        (
            const alphatWallFunctionFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new alphatWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
