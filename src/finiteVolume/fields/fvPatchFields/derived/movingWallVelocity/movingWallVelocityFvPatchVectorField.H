/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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
    Foam::movingWallVelocityFvPatchVectorField

Description
    This boundary condition provides a no-slip velocity condition for ridged
    moving walls or flexible moving walls on which the mesh vertices move
    with the surface.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            movingWallVelocity;
        value           uniform (0 0 0);    // Initial value
    }
    \endverbatim

See also
    Foam::fixedValueFvPatchVectorField

SourceFiles
    movingWallVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef movingWallVelocityFvPatchVectorField_H
#define movingWallVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class movingWallVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class movingWallVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{

public:

    //- Runtime type information
    TypeName("movingWallVelocity");


    // Constructors

        //- Construct from patch, internal field and dictionary
        movingWallVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given movingWallVelocityFvPatchVectorField
        //  onto a new patch
        movingWallVelocityFvPatchVectorField
        (
            const movingWallVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        movingWallVelocityFvPatchVectorField
        (
            const movingWallVelocityFvPatchVectorField&
        ) = delete;

        //- Copy constructor setting internal field reference
        movingWallVelocityFvPatchVectorField
        (
            const movingWallVelocityFvPatchVectorField&,
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
                new movingWallVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member Operators

        //- Inherit assignment
        using fixedValueFvPatchVectorField::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
