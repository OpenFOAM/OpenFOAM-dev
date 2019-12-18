/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::fixedPressureCompressibleDensityFvPatchScalarField

Description
    This boundary condition calculates a (liquid) compressible density as a
    function of pressure and fluid properties:

        \f[
            \rho = \rho_{l,sat} + \psi_l*(p - p_{sat})
        \f]

    where

    \vartable
        \rho    | density [kg/m^3]
        \rho_{l,sat} | saturation liquid density [kg/m^3]
        \psi_l  | liquid compressibility
        p       | pressure [Pa]
        p_{sat} | saturation pressure [Pa]
    \endvartable

    The variables \c rholSat, \c pSat and \c psil are retrieved from the
    \c thermodynamicProperties dictionary.

Usage
    \table
        Property     | Description             | Required    | Default value
        p            | pressure field name     | no          | p
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type        fixedPressureCompressibleDensity;
        p           p;
        value       uniform 1;
    }
    \endverbatim

See also
    Foam::fixedValueFvPatchField

SourceFiles
    fixedPressureCompressibleDensityFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedPressureCompressibleDensityFvPatchScalarField_H
#define fixedPressureCompressibleDensityFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class fixedPressureCompressibleDensityFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class fixedPressureCompressibleDensityFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private Data

        //- Name of static pressure field
        word pName_;


public:

    //- Runtime type information
    TypeName("fixedPressureCompressibleDensity");


    // Constructors

        //- Construct from patch and internal field
        fixedPressureCompressibleDensityFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedPressureCompressibleDensityFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  fixedPressureCompressibleDensityFvPatchScalarField
        //  onto a new patch
        fixedPressureCompressibleDensityFvPatchScalarField
        (
            const fixedPressureCompressibleDensityFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        fixedPressureCompressibleDensityFvPatchScalarField
        (
            const fixedPressureCompressibleDensityFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new fixedPressureCompressibleDensityFvPatchScalarField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        fixedPressureCompressibleDensityFvPatchScalarField
        (
            const fixedPressureCompressibleDensityFvPatchScalarField&,
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
                new fixedPressureCompressibleDensityFvPatchScalarField
                (
                    *this,
                    iF
                )
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
