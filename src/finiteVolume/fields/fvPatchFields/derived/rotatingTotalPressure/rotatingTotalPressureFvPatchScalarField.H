/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::rotatingTotalPressureFvPatchScalarField

Description
    This boundary condition provides a total pressure condition for patches
    in a rotating frame.

Usage
    \table
        Property     | Description             | Required    | Default value
        U            | velocity field name     | no          | U
        phi          | flux field name         | no          | phi
        rho          | density field name      | no          | none
        psi          | compressibility field name | no       | none
        gamma        | ratio of specific heats (Cp/Cv) | yes |
        p0           | static pressure reference | yes       |
        origin       | origin of rotation in Cartesian co-ordinates | yes|
        axis         | axis of rotation        | yes         |
        omega        | angular velocity of the frame [rad/s] | no |
        rpm          | angular velocity of the frame [rpm]   | no |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            rotatingTotalPressure;
        U               U;
        phi             phi;
        rho             rho;
        psi             psi;
        gamma           1.4;
        p0              uniform 1e5;
        origin          (0 0 0);
        axis            (0 0 1);
        rpm             100;
    }
    \endverbatim

    The \c omega or \c rpm entries are Function1 of time, see Foam::Function1s.

See also
    Foam::totalPressureFvPatchScalarField
    Foam::Function1s

SourceFiles
    rotatingTotalPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef rotatingTotalPressureFvPatchScalarField_H
#define rotatingTotalPressureFvPatchScalarField_H

#include "totalPressureFvPatchScalarField.H"
#include "omega.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class rotatingTotalPressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class rotatingTotalPressureFvPatchScalarField
:
    public totalPressureFvPatchScalarField
{
    // Private Data

        //- Origin of the rotation
        vector origin_;

        //- Axis of the rotation
        vector axis_;

        //- Angular velocity of the frame (rad/sec)
        Function1s::omega omega_;


public:

    //- Runtime type information
    TypeName("rotatingTotalPressure");


    // Constructors

        //- Construct from patch and internal field
        rotatingTotalPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        rotatingTotalPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given rotatingTotalPressureFvPatchScalarField
        //  onto a new patch
        rotatingTotalPressureFvPatchScalarField
        (
            const rotatingTotalPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        rotatingTotalPressureFvPatchScalarField
        (
            const rotatingTotalPressureFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        rotatingTotalPressureFvPatchScalarField
        (
            const rotatingTotalPressureFvPatchScalarField&,
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
                new rotatingTotalPressureFvPatchScalarField(*this, iF)
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
