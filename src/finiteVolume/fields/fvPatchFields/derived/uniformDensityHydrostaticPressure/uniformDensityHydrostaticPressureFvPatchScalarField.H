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
    Foam::uniformDensityHydrostaticPressureFvPatchScalarField

Description
    This boundary condition provides a hydrostatic pressure condition,
    calculated as:

        \f[
            p_{hyd} = p_{ref} + \rho_{ref} g (x - x_{ref})
        \f]

    where
    \vartable
        p_{hyd}    | Hydrostatic pressure [Pa]
        p_{ref}    | Reference pressure [Pa]
        x_{ref}    | Reference point in Cartesian co-ordinates
        \rho_{ref} | Density (assumed uniform)
        g          | Acceleration due to gravity [m/s^2]
    \endtable

Usage
    \table
        Property     | Description             | Required    | Default value
        rhoRef       | Uniform density [kg/m^3] | yes         |
        pRef         | Reference pressure [Pa] | yes         |
        pRefPoint    | Reference pressure location | no      | hRef
        value        | Initial value           | no          | pRef
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            uniformDensityHydrostaticPressure;
        rhoRef          1000;
        pRef            1e5;
        pRefPoint       (0 0 0);
    }
    \endverbatim

SourceFiles
    uniformDensityHydrostaticPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef uniformDensityHydrostaticPressureFvPatchScalarField_H
#define uniformDensityHydrostaticPressureFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
     Class uniformDensityHydrostaticPressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class uniformDensityHydrostaticPressureFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private Data

        //- Constant density in the far-field
        scalar rho_;

        //- Reference pressure
        scalar pRef_;

        //- True if the reference pressure location is specified
        bool pRefPointSpecified_;

        //- Optional reference pressure location
        vector pRefPoint_;


public:

    //- Runtime type information
    TypeName("uniformDensityHydrostaticPressure");


    // Constructors

        //- Construct from patch, internal field and dictionary
        uniformDensityHydrostaticPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  uniformDensityHydrostaticPressureFvPatchScalarField onto a new patch
        uniformDensityHydrostaticPressureFvPatchScalarField
        (
            const uniformDensityHydrostaticPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        uniformDensityHydrostaticPressureFvPatchScalarField
        (
            const uniformDensityHydrostaticPressureFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        uniformDensityHydrostaticPressureFvPatchScalarField
        (
            const uniformDensityHydrostaticPressureFvPatchScalarField&,
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
                new uniformDensityHydrostaticPressureFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        // Access

            //- Return the constant density in the far-field
            scalar rho() const
            {
                return rho_;
            }

            //- Return reference to the constant density in the far-field
            //  to allow adjustment
            scalar& rho()
            {
                return rho_;
            }

            //- Return the reference pressure
            scalar pRef() const
            {
                return pRef_;
            }

            //- Return reference to the reference pressure to allow adjustment
            scalar& pRef()
            {
                return pRef_;
            }

            //- Return the pressure reference location
            const vector& pRefPoint() const
            {
                return pRefPoint_;
            }

            //- Return reference to the pressure reference location
            //  to allow adjustment
            vector& pRefPoint()
            {
                return pRefPoint_;
            }


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
