/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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
    Foam::compressible::alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField

Group
    grpCmpWallFunctions

Description
    A simple alphatPhaseChangeWallFunctionFvPatchScalarField with a fixed
    volumetric phase-change mass flux.

SeeAlso
    Foam::compressible::alphatPhaseChangeWallFunctionFvPatchScalarField

SourceFiles
    alphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleAlphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField_H
#define compressibleAlphatFixedDmdtWallBoilingWallFunctionFvPatchScalarField_H

#include "alphatPhaseChangeWallFunctionFvPatchScalarField.H"

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
    public alphatPhaseChangeWallFunctionFvPatchScalarField
{
    // Private data

        //- Turbulent Prandtl number
        scalar Prt_;

        //- Cmu coefficient
        scalar Cmu_;

        //- Von Karman constant
        scalar kappa_;

        //- E coefficient
        scalar E_;

        //- E coefficient
        scalar fixedDmdt_;


        // Solution parameters

            static scalar maxExp_;
            static scalar tolerance_;
            static label maxIters_;


    // Private Member Functions

        //- Check the type of the patch
        void checkType();

        //- 'P' function
        tmp<scalarField> Psmooth(const scalarField& Prat) const;

        //- Calculate y+ at the edge of the thermal laminar sublayer
        tmp<scalarField> yPlusTherm
        (
            const scalarField& P,
            const scalarField& Prat
        ) const;


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

        //- Return the rate of phase-change
        virtual const scalarField& dmdt() const
        {
            return dmdt_;
        }


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
