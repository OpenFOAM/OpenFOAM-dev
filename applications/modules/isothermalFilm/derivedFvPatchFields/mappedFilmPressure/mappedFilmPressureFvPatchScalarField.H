/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::mappedFilmPressureFvPatchScalarField

Description
    Film pressure boundary condition which maps the neighbouring fluid patch
    pressure to both the surface patch and internal film pressure field.

SourceFiles
    mappedFilmPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef mappedFilmPressureFvPatchScalarField_H
#define mappedFilmPressureFvPatchScalarField_H

#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class mappedFilmPressureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class mappedFilmPressureFvPatchScalarField
:
    public zeroGradientFvPatchField<scalar>
{

public:

    //- Runtime type information
    TypeName("mappedFilmPressure");


    // Constructors

        //- Construct from patch, internal field and dictionary
        mappedFilmPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given mappedFilmPressureFvPatchScalarField
        //  onto a new patch
        mappedFilmPressureFvPatchScalarField
        (
            const mappedFilmPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        mappedFilmPressureFvPatchScalarField
        (
            const mappedFilmPressureFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        mappedFilmPressureFvPatchScalarField
        (
            const mappedFilmPressureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<scalar>>
            (
                new mappedFilmPressureFvPatchScalarField(*this, iF)
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
