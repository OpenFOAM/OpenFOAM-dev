/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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
    Foam::movingWallSlipVelocityFvPatchVectorField

Description
    This boundary condition provides a slip velocity condition for cases with
    moving walls.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            movingWallSlipVelocity;
        value           uniform (0 0 0);    // Initial value
    }
    \endverbatim

SourceFiles
    movingWallSlipVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef movingWallSlipVelocityFvPatchVectorField_H
#define movingWallSlipVelocityFvPatchVectorField_H

#include "directionMixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class movingWallSlipVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class movingWallSlipVelocityFvPatchVectorField
:
    public directionMixedFvPatchVectorField
{

public:

    //- Runtime type information
    TypeName("movingWallSlipVelocity");


    // Constructors

        //- Construct from patch and internal field
        movingWallSlipVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        movingWallSlipVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given movingWallSlipVelocityFvPatchVectorField
        //  onto a new patch
        movingWallSlipVelocityFvPatchVectorField
        (
            const movingWallSlipVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        movingWallSlipVelocityFvPatchVectorField
        (
            const movingWallSlipVelocityFvPatchVectorField&
        ) = delete;

        //- Copy constructor setting internal field reference
        movingWallSlipVelocityFvPatchVectorField
        (
            const movingWallSlipVelocityFvPatchVectorField&,
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
                new movingWallSlipVelocityFvPatchVectorField(*this, iF)
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
