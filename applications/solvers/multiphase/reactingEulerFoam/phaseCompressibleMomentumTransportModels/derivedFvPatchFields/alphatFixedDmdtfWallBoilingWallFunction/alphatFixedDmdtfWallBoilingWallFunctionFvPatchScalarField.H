/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2020 OpenFOAM Foundation
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
    Foam::compressible::
        alphatFixedDmdtfWallBoilingWallFunctionFvPatchScalarField

Description
    A simple alphatPhaseChangeWallFunctionFvPatchScalarField with
    a fixed volumetric phase-change mass flux.

See also
    Foam::compressible::alphatPhaseChangeWallFunctionFvPatchScalarField

SourceFiles
    alphatFixedDmdtfWallBoilingWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef alphatFixedDmdtfWallBoilingWallFunctionFvPatchScalarField_H
#define alphatFixedDmdtfWallBoilingWallFunctionFvPatchScalarField_H

#include "alphatPhaseChangeWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
   Class alphatFixedDmdtfWallBoilingWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class alphatFixedDmdtfWallBoilingWallFunctionFvPatchScalarField
:
    public alphatPhaseChangeWallFunctionFvPatchScalarField
{
    // Private Data

        //- Volumetric phase-change mass flux in near wall cells
        scalar fixedDmdtf_;

        //- Latent heat
        scalar L_;


public:

    //- Runtime type information
    TypeName("compressible::alphatFixedDmdtfWallBoilingWallFunction");


    // Constructors

        //- Construct from patch and internal field
        alphatFixedDmdtfWallBoilingWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        alphatFixedDmdtfWallBoilingWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  alphatFixedDmdtfWallBoilingWallFunctionFvPatchScalarField
        //  onto a new patch
        alphatFixedDmdtfWallBoilingWallFunctionFvPatchScalarField
        (
            const alphatFixedDmdtfWallBoilingWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        alphatFixedDmdtfWallBoilingWallFunctionFvPatchScalarField
        (
            const alphatFixedDmdtfWallBoilingWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new alphatFixedDmdtfWallBoilingWallFunctionFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Copy constructor setting internal field reference
        alphatFixedDmdtfWallBoilingWallFunctionFvPatchScalarField
        (
            const alphatFixedDmdtfWallBoilingWallFunctionFvPatchScalarField&,
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
                new alphatFixedDmdtfWallBoilingWallFunctionFvPatchScalarField
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


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
