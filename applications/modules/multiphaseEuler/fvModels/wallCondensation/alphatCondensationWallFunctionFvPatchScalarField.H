/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::alphatCondensationWallFunctionFvPatchScalarField

Description
    This boundary condition provides a thermal wall function or the turbulent
    thermal diffusivity for walls on which a wall condensation fvModel is in
    operation. The controls for the condensation model are not specified in this
    condition; they are provided to the wall condensation fvModel.
    The presence of this condition is used, however, to determine on which walls
    condensation takes place.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchType>
    {
        type            alphatCondensationWallFunction;
        libs            ("libmultiphaseEulerFvModels.so");
        value           uniform 0;
    }
    \endverbatim

See also
    Foam::fv::wallCondensation

SourceFiles
    alphatCondensationWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef alphatCondensationWallFunctionFvPatchScalarField_H
#define alphatCondensationWallFunctionFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "wallCondensation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class alphatCondensationWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class alphatCondensationWallFunctionFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
private:

    // Private Member Data

        //- The name of the wall condensation model
        const word modelName_;


    // Private Member Functions

        //- Access the wall condensation model
        const fv::wallCondensation& model() const;


public:

    //- Runtime type information
    TypeName("alphatCondensationWallFunction");


    // Constructors

        //- Construct from patch, internal field and dictionary
        alphatCondensationWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  alphatCondensationWallFunctionFvPatchScalarField
        //  onto a new patch
        alphatCondensationWallFunctionFvPatchScalarField
        (
            const alphatCondensationWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        alphatCondensationWallFunctionFvPatchScalarField
        (
            const alphatCondensationWallFunctionFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        alphatCondensationWallFunctionFvPatchScalarField
        (
            const alphatCondensationWallFunctionFvPatchScalarField&,
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
                new alphatCondensationWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            virtual void write(Ostream&) const;


    // Member Operators

        //- Inherit assignment
        using fixedValueFvPatchScalarField::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
