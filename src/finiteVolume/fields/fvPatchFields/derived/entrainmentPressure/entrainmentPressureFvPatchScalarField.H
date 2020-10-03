/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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
    Foam::entrainmentPressureFvPatchScalarField

Description
    This boundary condition provides an entrainment pressure condition. The
    condition is the same as totalPressure condition except that where the
    total pressure uses the square of the velocity magnitude, this uses the
    signed square of the normal velocity. In the incompressible case, that
    means the formulation is as follows:

        \f[
            p_p = p_0 + 0.5 Un |Un|
        \f]
        where
        \vartable
            p_p     | pressure at patch [m^2/s^2]
            p_0     | reference pressure [m^2/s^2]
            U       | patch normal velocity [m/s]
        \endvartable

    This condition is particularly suitable for boundaries on which the flow
    direction changes frequently and/or when the inlet/outlet direction is not
    well defined.

Usage
    \table
        Property     | Description                | Required | Default value
        phi          | Flux field name            | no       | phi
        rho          | Density field name         | no       | rho
        psi          | Compressibility field name | no       | none
        gamma        | (Cp/Cv)                    | no       | 1
        p0           | Reference pressure         | yes      |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            entrainmentPressure;
        p0              uniform 1e5;
    }
    \endverbatim

See also
    Foam::dynamicPressureFvPatchScalarField
    Foam::totalPressureFvPatchScalarField
    Foam::fixedValueFvPatchField

SourceFiles
    entrainmentPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef entrainmentPressureFvPatchScalarField_H
#define entrainmentPressureFvPatchScalarField_H

#include "dynamicPressureFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class entrainmentPressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class entrainmentPressureFvPatchScalarField
:
    public dynamicPressureFvPatchScalarField
{
    // Private Data

        //- Name of the flux field
        const word phiName_;


public:

    //- Runtime type information
    TypeName("entrainmentPressure");


    // Constructors

        //- Construct from patch and internal field
        entrainmentPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        entrainmentPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given entrainmentPressureFvPatchScalarField
        //  onto a new patch
        entrainmentPressureFvPatchScalarField
        (
            const entrainmentPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        entrainmentPressureFvPatchScalarField
        (
            const entrainmentPressureFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        entrainmentPressureFvPatchScalarField
        (
            const entrainmentPressureFvPatchScalarField&,
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
                new entrainmentPressureFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

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
