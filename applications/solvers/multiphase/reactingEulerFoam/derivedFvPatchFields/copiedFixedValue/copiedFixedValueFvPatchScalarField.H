/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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
    Foam::copiedFixedValueFvPatchScalarField

Description
    Copies the boundary values from a user specified field.

See also
    Foam::fixedValueFvPatchField

SourceFiles
    copiedFixedValueFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef copiedFixedValueFvPatchScalarField_H
#define copiedFixedValueFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class copiedFixedValueFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class copiedFixedValueFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
protected:

    // Protected data

        word sourceFieldName_;

public:

    //- Runtime type information
    TypeName("copiedFixedValue");


    // Constructors

        //- Construct from patch and internal field
        copiedFixedValueFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        copiedFixedValueFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  copiedFixedValueFvPatchScalarField
        //  onto a new patch
        copiedFixedValueFvPatchScalarField
        (
            const copiedFixedValueFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        copiedFixedValueFvPatchScalarField
        (
            const copiedFixedValueFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        copiedFixedValueFvPatchScalarField
        (
            const copiedFixedValueFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
