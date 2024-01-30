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
    Foam::movingMappedWallVelocityFvPatchVectorField

Description
    This boundary condition provides a no-slip velocity condition for mapped
    walls. The wall velocity is taken to be the mesh velocity of the
    neighbouring region.

    This will typically be used in CHT simulations in order to apply the
    mesh motion of a solid region to the boundary of the adjacent fluid
    region.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            movingMappedWallVelocity;
        value           uniform (0 0 0);    // Initial value
    }
    \endverbatim

See also
    Foam::fixedValueFvPatchVectorField

SourceFiles
    movingMappedWallVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef movingMappedWallVelocityFvPatchVectorField_H
#define movingMappedWallVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class movingMappedWallVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class movingMappedWallVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{

public:

    //- Runtime type information
    TypeName("movingMappedWallVelocity");


    // Constructors

        //- Construct from patch, internal field and dictionary
        movingMappedWallVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  movingMappedWallVelocityFvPatchVectorField onto a new patch
        movingMappedWallVelocityFvPatchVectorField
        (
            const movingMappedWallVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        movingMappedWallVelocityFvPatchVectorField
        (
            const movingMappedWallVelocityFvPatchVectorField&
        ) = delete;

        //- Copy constructor setting internal field reference
        movingMappedWallVelocityFvPatchVectorField
        (
            const movingMappedWallVelocityFvPatchVectorField&,
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
                new movingMappedWallVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member Functions

        // Mapping functions

            //- Map the given fvPatchField onto this fvPatchField
            virtual void map(const fvPatchVectorField&, const fieldMapper&);


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
