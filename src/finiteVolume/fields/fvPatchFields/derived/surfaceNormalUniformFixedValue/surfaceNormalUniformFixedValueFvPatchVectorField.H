/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2020 OpenFOAM Foundation
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
    Foam::surfaceNormalUniformFixedValueFvPatchVectorField

Description
    Surface-normal fixed value vector boundary condition

    in which the vector field is obtained by scaling the patch normals
    by the evaluated scalar Function1 provided.

Usage
    \table
        Property     | Description             | Required    | Default value
        uniformValue | uniform value           | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            surfaceNormalUniformFixedValue;
        uniformValue    constant -0.2;         // 0.2 INTO the domain
    }
    \endverbatim

Note
    Sign conventions:
    - the value is positive for outward-pointing vectors

See also
    Foam::fixedValueFvPatchField

SourceFiles
    surfaceNormalUniformFixedValueFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceNormalUniformFixedValueFvPatchVectorField_H
#define surfaceNormalUniformFixedValueFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class surfaceNormalUniformFixedValueFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class surfaceNormalUniformFixedValueFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private Data

        autoPtr<Function1<scalar>> uniformValue_;


public:

    //- Runtime type information
    TypeName("surfaceNormalUniformFixedValue");


    // Constructors

        //- Construct from patch and internal field
        surfaceNormalUniformFixedValueFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        surfaceNormalUniformFixedValueFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  surfaceNormalUniformFixedValueFvPatchVectorField
        //  onto a new patch
        surfaceNormalUniformFixedValueFvPatchVectorField
        (
            const surfaceNormalUniformFixedValueFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        surfaceNormalUniformFixedValueFvPatchVectorField
        (
            const surfaceNormalUniformFixedValueFvPatchVectorField&
        ) = delete;

        //- Copy constructor setting internal field reference
        surfaceNormalUniformFixedValueFvPatchVectorField
        (
            const surfaceNormalUniformFixedValueFvPatchVectorField&,
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
                new surfaceNormalUniformFixedValueFvPatchVectorField
                (
                    *this,
                    iF
                )
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
