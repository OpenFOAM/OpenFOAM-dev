/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2024 OpenFOAM Foundation
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
    Foam::swirlInletVelocityFvPatchVectorField

Description
    Velocity inlet boundary condition creating axial, radial and tangential
    velocity fields specified by functions of time and radius or by a given
    angular speed.

Usage
    \table
        Property     | Description             | Required    | Default value
        origin       | Origin of rotation      | no          | patch centre
        axis         | Axis of rotation        | no          | patch normal
        axialVelocity| Axial velocity [m/s]    | yes         |
        radialVelocity | Radial velocity [m/s] | yes         |
        omega        | Angular velocity of the frame | no |
        tangentialVelocity | Tangential velocity [m/s] | no  |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            swirlInletVelocity;
        origin          (0 0 0);
        axis            (0 0 1);
        axialVelocity   constant 30;
        radialVelocity  constant 10;
        tangentialVelocity  constant 100;
    }
    \endverbatim

    The \c axialVelocity, \c radialVelocity and \c tangentialVelocity
    entries are Function2 of time and radius, see see Foam::Function2s. The \c
    omega is a Function1 of time, see Foam::Function1s.

See also
    Foam::fixedValueFvPatchField
    Foam::Function1s
    Foam::Function2s

SourceFiles
    swirlInletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef swirlInletVelocityFvPatchVectorField_H
#define swirlInletVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "Function2.H"
#include "omega1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
            Class swirlInletVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class swirlInletVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private Data

        //- Origin of the rotation
        const vector origin_;

        //- Axis of the rotation
        const vector axis_;

        //- Axial velocity
        autoPtr<Function2<scalar>> axialVelocity_;

        //- Radial velocity
        autoPtr<Function2<scalar>> radialVelocity_;

        //- Angular velocity of the frame
        autoPtr<Function1s::omega> omega_;

        //- Tangential velocity
        autoPtr<Function2<scalar>> tangentialVelocity_;


public:

   //- Runtime type information
   TypeName("swirlInletVelocity");


   // Constructors

        //- Construct from patch, internal field and dictionary
        swirlInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  flowRateInletVelocityFvPatchVectorField
        //  onto a new patch
        swirlInletVelocityFvPatchVectorField
        (
            const swirlInletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        swirlInletVelocityFvPatchVectorField
        (
            const swirlInletVelocityFvPatchVectorField&
        ) = delete;

        //- Copy constructor setting internal field reference
        swirlInletVelocityFvPatchVectorField
        (
            const swirlInletVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new swirlInletVelocityFvPatchVectorField(*this, iF)
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
