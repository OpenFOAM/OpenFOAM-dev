/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 OpenFOAM Foundation
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
        alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField

Group
    grpCmpWallFunctions

Description
    A simple alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField with
    a fixed volumetric phase-change mass flux.

See also
    Foam::compressible::
        alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField

SourceFiles
    alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleAlphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField_H
#define compressibleAlphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField_H

#include "alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
   Class alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField
:
    public alphatPhaseChangeJayatillekeWallFunctionFvPatchScalarField
{
    // Private data

        //- dmdt relaxationFactor
        scalar relax_;

        //- Volumetric phase-change mass flux in near wall cells
        scalar fixedDmdt_;

        //- Latent heat
        scalar L_;


public:

    //- Runtime type information
    TypeName("compressible::alphatFixedDmdtWallBoilingWallFunction");


    // Constructors

        //- Construct from patch and internal field
        alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField
        //  onto a new patch
        alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField
        (
            const alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField
        (
            const alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField
        (
            const alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField&,
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
                new alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

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
