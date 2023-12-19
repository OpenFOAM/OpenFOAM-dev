/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2023 OpenFOAM Foundation
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
    This is a boundary condition for pressure suitable for boundaries at which
    the flow direction is uncertain but both inflow and/or outflow can occur.
    The condition switches based on the direction of flow. For outflow, the
    patch pressure is simply set to the reference pressure. For inflow, the
    patch pressure is reduced by the dynamic pressure calculated using the
    inflow velocity normal to the patch.

    For incompressible flow, with pressure expressed in kinematic units, the
    pressure with inflow is calculated as:

        \f[
            p_p = p_0 - 0.5 |Un|^2
        \f]
        where
        \vartable
            p_p     | pressure at patch [m^2/s^2]
            p_0     | reference pressure [m^2/s^2]
            Un      | patch normal velocity [m/s]
        \endvartable

    The condition is similar to the totalPressure condition and includes
    further calculations of pressure for compressible flow under subsonic,
    transonic and supersonic regimes.

    Like totalPressure, entrainmentPressure is usually applied in conjunction
    with the pressureInletOutletVelocity condition for velocity. The
    entrainmentPressure condition calculates the pressure using a normal
    velocity, based on the fluxes at the boundary. The totalPressure condition
    calculates the pressure using the inflow velocity itself.

    With the pressureInletOutletVelocity condition, using the default
    tangentialVelocity of zero, any inflow velocity is constrained to be
    normal to the boundary, making the two pressure boundary conditions
    very similar. The entrainmentPressure condition is, however, more robust
    because it couples with the flux rather than velocity. The condition
    performs particularly well at boundaries where the flow direction tends
    to switch frequently between in and out of the domain.

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
            const fieldMapper&
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
