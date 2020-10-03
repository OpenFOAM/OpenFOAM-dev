/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2020 OpenFOAM Foundation
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
    This boundary condition describes an inlet vector boundary condition in
    swirl co-ordinates given a central axis, central point, axial, radial and
    tangential velocity profiles.

Usage
    \table
        Property       | Description             | Required    | Default value
        axis           | Axis of rotation        | yes         |
        origin         | Origin of rotation      | yes         |
        axialVelocity  | Axial velocity profile [m/s] | yes    |
        radialVelocity | Radial velocity profile [m/s] | yes  |
        tangentialVelocity | tangential velocity profile [m/s] | yes  |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            swirlInletVelocity;
        axis            (0 0 1);
        origin          (0 0 0);
        axialVelocity   constant 30;
        radialVelocity  constant 10;
        tangentialVelocity  constant 100;
    }
    \endverbatim

Note
    The \c axialVelocity, \c radialVelocity and \c tangentialVelocity entries
    are Function1 types, able to describe time varying functions.  The example
    above gives the usage for supplying constant values.

See also
    Foam::fixedValueFvPatchField
    Foam::Function1s

SourceFiles
    swirlInletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef swirlInletVelocityFvPatchVectorField_H
#define swirlInletVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "Function1.H"

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
        autoPtr<Function1<scalar>> axialVelocity_;

        //- Radial velocity
        autoPtr<Function1<scalar>> radialVelocity_;

        //- Tangential velocity
        autoPtr<Function1<scalar>> tangentialVelocity_;


public:

   //- Runtime type information
   TypeName("swirlInletVelocity");


   // Constructors

        //- Construct from patch and internal field
        swirlInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

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
            const fvPatchFieldMapper&
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
