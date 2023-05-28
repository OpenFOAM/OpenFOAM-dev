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
    Foam::totalPressureFvPatchScalarField

Description
    Inflow, outflow and entrainment pressure boundary condition based on a
    constant total pressure assumption.

    For outflow the patch pressure is set to the external static pressure.

    For inflow the patch pressure is evaluated from the patch velocity and the
    external total pressure obtained from the external static pressure \c p_0
    and external velocity \c U_0 which is looked-up from the the optional \c
    tangentialVelocity entry in the \c pressureInletOutletVelocity velocity
    boundary condition for the patch if that boundary condition is used,
    otherwise \c U_0 is assumed zero and the external total pressure is equal to
    the external static pressure.

    The patch pressure is evaluated from the external conditions using one of
    the following expressions depending on the flow conditions and
    specification of compressibility:

    1. incompressible subsonic:
        \f[
            p_p = p_0 + 0.5 |U_0|^2 - 0.5 |U|^2
        \f]
        where
        \vartable
            p_p     | pressure at patch [m^2/s^2]
            p_0     | external static pressure [m^2/s^2]
            U       | velocity [m/s]
            U_0     | external velocity [m/s]
        \endvartable

    2. compressible subsonic:
        \f[
            p_p = p_0 + \rho (0.5 |U_0|^2 - 0.5 |U|^2)
        \f]
        where
        \vartable
            p_p     | pressure at patch [Pa]
            p_0     | external static pressure [Pa]
            \rho    | density [kg/m^3]
            U       | velocity [m/s]
            U_0     | external velocity [m/s]
        \endvartable

    3. compressible transonic (\f$\gamma = 1\f$):
        \f[
            p_p = \frac{p_0}{1 + \psi (0.5 |U|^2 - 0.5 |U_0|^2)}
        \f]
        where
        \vartable
            p_p     | pressure at patch [Pa]
            p_0     | external static pressure [Pa]
            \psi    | compressibility [m^2/s^2]
            \rho    | density [kg/m^3]
            U       | velocity [m/s]
            U_0     | external velocity [m/s]
        \endvartable

    4. compressible supersonic (\f$\gamma > 1\f$):
        \f[
            p_p = \frac{p_0}
                  {(1 +  G \psi (0.5 |U|^2 - 0.5 |U_0|^2))^{\frac{1}{G}}}
        \f]
        where
        \vartable
            p_p     | pressure at patch [Pa]
            p_0     | external static pressure [Pa]
            \psi    | compressibility [m^2/s^2]
            \rho    | density [kg/m^3]
            G       | coefficient given by \f$\frac{\gamma}{1-\gamma}\f$ []
            \gamma  | ratio of specific heats (Cp/Cv) []
            U       | velocity [m/s]
            U_0     | external velocity [m/s]
        \endvartable

    The modes of operation are set by the dimensions of the pressure field
    to which this boundary condition is applied, the \c psi entry and the value
    of \c gamma:
    \table
        Mode                    | dimensions | psi   | gamma
        incompressible subsonic | p/rho      |       |
        compressible subsonic   | p          | none  |
        compressible transonic  | p          | psi   | 1
        compressible supersonic | p          | psi   | > 1
    \endtable

Usage
    \table
        Property     | Description                | Required | Default value
        U            | Velocity field name        | no       | U
        phi          | Flux field name            | no       | phi
        rho          | Density field name         | no       | rho
        psi          | Compressibility field name | no       | none
        gamma        | (Cp/Cv)                    | no       | 1
        p0           | External pressure          | yes      |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            totalPressure;
        p0              uniform 1e5;
    }
    \endverbatim

See also
    Foam::dynamicPressureFvPatchScalarField
    Foam::fixedValueFvPatchField
    Foam::pressureInletOutletVelocityFvPatchVectorField

SourceFiles
    totalPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef totalPressureFvPatchScalarField_H
#define totalPressureFvPatchScalarField_H

#include "dynamicPressureFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class totalPressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class totalPressureFvPatchScalarField
:
    public dynamicPressureFvPatchScalarField
{
protected:

    // Protected Data

        //- Name of the velocity field
        const word UName_;

        //- Name of the flux field
        const word phiName_;


public:

    //- Runtime type information
    TypeName("totalPressure");


    // Constructors

        //- Construct from patch, internal field and dictionary
        totalPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given totalPressureFvPatchScalarField
        //  onto a new patch
        totalPressureFvPatchScalarField
        (
            const totalPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        totalPressureFvPatchScalarField
        (
            const totalPressureFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        totalPressureFvPatchScalarField
        (
            const totalPressureFvPatchScalarField&,
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
                new totalPressureFvPatchScalarField(*this, iF)
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
