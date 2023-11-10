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
    Foam::phaseHydrostaticPressureFvPatchScalarField

Description
    This boundary condition provides a phase-based hydrostatic pressure
    condition, calculated as:

        \f[
            p_{hyd} = p_{ref} + \rho g (x - x_{ref})
        \f]

    where
    \vartable
        p_{hyd} | hyrostatic pressure [Pa]
        p_{ref} | reference pressure [Pa]
        x_{ref} | reference point in Cartesian co-ordinates
        \rho    | density (assumed uniform)
        g       | acceleration due to gravity [m/s^2]
    \endtable

    The values are assigned according to the phase-fraction field:
    - 1: apply \$fp_{hyd}\$f
    - 0: apply a zero-gradient condition

Usage
    \table
        Property      | Description                 | Required | Default value
        phaseFraction | phase-fraction field name   | no       | alpha
        rho           | density field name          | no       | rho
        pRefValue     | reference pressure [Pa]     | yes      |
        pRefPoint     | reference pressure location | yes      |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            phaseHydrostaticPressure;
        phaseFraction   alpha1;
        rho             rho;
        pRefValue       1e5;
        pRefPoint       (0 0 0);
        value           uniform 0; // optional initial value
    }
    \endverbatim

See also
    Foam::mixedFvPatchScalarField

SourceFiles
    phaseHydrostaticPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef phaseHydrostaticPressureFvPatchScalarField_H
#define phaseHydrostaticPressureFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class phaseHydrostaticPressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class phaseHydrostaticPressureFvPatchScalarField
:
    public mixedFvPatchScalarField
{

protected:

    // Protected data

        //- Name of phase-fraction field
        word phaseFraction_;

        //- Constant density in the far-field
        scalar rho_;

        //- Reference pressure
        scalar pRefValue_;

        //- Reference pressure location
        vector pRefPoint_;


public:

    //- Runtime type information
    TypeName("phaseHydrostaticPressure");


    // Constructors

        //- Construct from patch, internal field and dictionary
        phaseHydrostaticPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  phaseHydrostaticPressureFvPatchScalarField onto a new patch
        phaseHydrostaticPressureFvPatchScalarField
        (
            const phaseHydrostaticPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        phaseHydrostaticPressureFvPatchScalarField
        (
            const phaseHydrostaticPressureFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        phaseHydrostaticPressureFvPatchScalarField
        (
            const phaseHydrostaticPressureFvPatchScalarField&,
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
                new phaseHydrostaticPressureFvPatchScalarField(*this, iF)
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

            //- Return the phaseFraction
            const word& phaseFraction() const
            {
                return phaseFraction_;
            }

            //- Return reference to the phaseFraction to allow adjustment
            word& phaseFraction()
            {
                return phaseFraction_;
            }

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


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member Operators

        virtual void operator=(const fvPatchScalarField& pvf);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
