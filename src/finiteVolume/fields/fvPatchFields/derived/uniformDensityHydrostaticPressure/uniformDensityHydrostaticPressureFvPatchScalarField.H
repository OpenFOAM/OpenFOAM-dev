/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides a hydrostatic pressure condition,
    calculated as:

        \f[
            p_{hyd} = p_{ref} + \rho g (x - x_{ref})
        \f]

    where
    \vartable
        p_{hyd} | Hydrostatic pressure [Pa]
        p_{ref} | Reference pressure [Pa]
        x_{ref} | Reference point in Cartesian co-ordinates
        \rho    | Density (assumed uniform)
        g       | Acceleration due to gravity [m/s2]
    \endtable

Usage
    \table
        Property     | Description             | Required    | Default value
        rho          | Uniform density [kg/m3] | yes         |
        pRefValue    | Reference pressure [Pa] | yes         |
        pRefPoint    | Reference pressure location | yes     |
        value        | Initial value           | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            uniformDensityHydrostaticPressure;
        rho             rho;
        pRefValue       1e5;
        pRefPoint       (0 0 0);
        value           uniform 0; // Required initial value
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
    // Private data

        //- Constant density in the far-field
        scalar rho_;

        //- Reference pressure
        scalar pRefValue_;

        //- Reference pressure location
        vector pRefPoint_;


public:

    //- Runtime type information
    TypeName("uniformDensityHydrostaticPressure");


    // Constructors

        //- Construct from patch and internal field
        uniformDensityHydrostaticPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

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
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        uniformDensityHydrostaticPressureFvPatchScalarField
        (
            const uniformDensityHydrostaticPressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new uniformDensityHydrostaticPressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
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


    // Member functions

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
            scalar pRefValue() const
            {
                return pRefValue_;
            }

            //- Return reference to the reference pressure to allow adjustment
            scalar& pRefValue()
            {
                return pRefValue_;
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
