/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::alphatBoilingWallFunctionFvPatchScalarField

Description
    This boundary condition provides a thermal wall function or the turbulent
    thermal diffusivity for walls on which a wall boiling fvModel is in
    operation. The controls for the boiling model are not specified in this
    condition; they are provided to the wall boiling fvModel. The presence of
    this condition is used, however, to determine on which walls boiling takes
    place.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchType>
    {
        type            alphatBoilingWallFunction;
        libs            ("libmultiphaseEulerFvModels.so");
        value           uniform 0;
    }
    \endverbatim

See also
    Foam::fv::wallBoiling

SourceFiles
    alphatBoilingWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef alphatBoilingWallFunctionFvPatchScalarField_H
#define alphatBoilingWallFunctionFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "wallBoiling.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class alphatBoilingWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class alphatBoilingWallFunctionFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
private:

    // Private Member Data

        //- The name of the wall boiling model
        const word modelName_;


    // Private Member Functions

        //- Access the wall boiling model
        const fv::wallBoiling& model() const;


public:

    //- Runtime type information
    TypeName("alphatBoilingWallFunction");


    // Constructors

        //- Construct from patch, internal field and dictionary
        alphatBoilingWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  alphatBoilingWallFunctionFvPatchScalarField
        //  onto a new patch
        alphatBoilingWallFunctionFvPatchScalarField
        (
            const alphatBoilingWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        alphatBoilingWallFunctionFvPatchScalarField
        (
            const alphatBoilingWallFunctionFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        alphatBoilingWallFunctionFvPatchScalarField
        (
            const alphatBoilingWallFunctionFvPatchScalarField&,
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
                new alphatBoilingWallFunctionFvPatchScalarField(*this, iF)
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

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
