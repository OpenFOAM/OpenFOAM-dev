/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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
    Foam::prghUniformDensityHydrostaticPressureFvPatchScalarField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides a hydrostatic pressure condition for p_rgh,
    calculated as:

        \f[
            p_{rgh} = p_{ref} - (\rho - \rho_0) g (h - h_{ref})
        \f]

    where
    \vartable
        p_{rgh}    | Pseudo hydrostatic pressure [Pa]
        p_{ref}    | Static pressure at hRef [Pa]
        h          | Height in the opposite direction to gravity
        h_{ref}    | Reference height in the opposite direction to gravity
        \rho       | Density field
        \rho_{ref} | Uniform reference density at boundary
        g          | Acceleration due to gravity [m/s^2]
    \endtable

Usage
    \table
        Property     | Description               | Required    | Default value
        pRef         | Reference static pressure | yes         |
        rhoRef       | Reference density         | yes         |
        rho          | Density field name        | no          | rho
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            prghUniformDensityHydrostaticPressure;
        rhoRef          1000;
        p               0;
        value           uniform 0; // optional initial value
    }
    \endverbatim

See also
    Foam::fixedValueFvPatchScalarField

SourceFiles
    prghUniformDensityHydrostaticPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef prghUniformDensityHydrostaticPressureFvPatchScalarField_H
#define prghUniformDensityHydrostaticPressureFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
   Class prghUniformDensityHydrostaticPressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class prghUniformDensityHydrostaticPressureFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{

protected:

    // Protected data

        //- Reference Static pressure
        scalar pRef_;

        //- Reference density
        scalar rhoRef_;

        //- Name of phase-fraction field
        word rhoName_;


public:

    //- Runtime type information
    TypeName("prghUniformDensityHydrostaticPressure");


    // Constructors

        //- Construct from patch and internal field
        prghUniformDensityHydrostaticPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        prghUniformDensityHydrostaticPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  prghUniformDensityHydrostaticPressureFvPatchScalarField
        //  onto a new patch
        prghUniformDensityHydrostaticPressureFvPatchScalarField
        (
            const prghUniformDensityHydrostaticPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        prghUniformDensityHydrostaticPressureFvPatchScalarField
        (
            const prghUniformDensityHydrostaticPressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField >
            (
                new prghUniformDensityHydrostaticPressureFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        prghUniformDensityHydrostaticPressureFvPatchScalarField
        (
            const prghUniformDensityHydrostaticPressureFvPatchScalarField&,
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
                new prghUniformDensityHydrostaticPressureFvPatchScalarField
                (
                    *this,
                    iF
                )
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
