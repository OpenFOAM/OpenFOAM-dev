/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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
    Foam::prghTotalHydrostaticPressureFvPatchScalarField

Description
    This boundary condition provides static pressure condition for p_rgh,
    calculated as:

        \f[
            p_rgh = ph_rgh - 0.5 \rho |U|^2
        \f]

    where
    \vartable
        p_rgh   | Pressure - \rho g.(h - hRef) [Pa]
        ph_rgh  | Hydrostatic pressure - \rho g.(h - hRef) [Pa]
        h       | Height in the opposite direction to gravity
        hRef    | Reference height in the opposite direction to gravity
        \rho    | Density
        g       | Acceleration due to gravity [m/s^2]
    \endtable

Usage
    \table
        Property     | Description             | Required    | Default value
        U            | Velocity field name     | no          | U
        phi          | Flux field name         | no          | phi
        rho          | Density field name      | no          | rho
        ph_rgh       | ph_rgh field name       | no          | ph_rgh
        value        | Patch face values       | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            prghTotalHydrostaticPressure;
        value           uniform 0;
    }
    \endverbatim

See also
    Foam::fixedValueFvPatchScalarField
    Foam::prghTotalPressureFvPatchScalarField

SourceFiles
    prghTotalHydrostaticPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef prghTotalHydrostaticPressureFvPatchScalarField_H
#define prghTotalHydrostaticPressureFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class prghTotalHydrostaticPressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class prghTotalHydrostaticPressureFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{

protected:

    // Protected data

        //- Name of the velocity field
        word UName_;

        //- Name of the flux transporting the field
        word phiName_;

        //- Name of density field
        word rhoName_;

        //- Name of hydrostatic pressure field
        word ph_rghName_;


public:

    //- Runtime type information
    TypeName("prghTotalHydrostaticPressure");


    // Constructors

        //- Construct from patch and internal field
        prghTotalHydrostaticPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        prghTotalHydrostaticPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  prghTotalHydrostaticPressureFvPatchScalarField onto a new patch
        prghTotalHydrostaticPressureFvPatchScalarField
        (
            const prghTotalHydrostaticPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        prghTotalHydrostaticPressureFvPatchScalarField
        (
            const prghTotalHydrostaticPressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField >
            (
                new prghTotalHydrostaticPressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        prghTotalHydrostaticPressureFvPatchScalarField
        (
            const prghTotalHydrostaticPressureFvPatchScalarField&,
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
                new prghTotalHydrostaticPressureFvPatchScalarField(*this, iF)
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
