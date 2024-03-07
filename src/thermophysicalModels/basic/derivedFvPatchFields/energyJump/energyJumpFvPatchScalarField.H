/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2024 OpenFOAM Foundation
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
    Foam::energyJumpFvPatchScalarField

Description
    This boundary condition provides an energy jump condition, using the \c
    cyclic condition as a base. The jump is calculated from the corresponding
    temperature jump. This condition is not applied directly by the user, but
    is constructed automatically by the thermodynamic model when the
    corresponding temperature boundary condition is a jump cyclic.

See also
    Foam::jumpCyclicFvPatchField

SourceFiles
    energyJumpFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef energyJumpFvPatchScalarField_H
#define energyJumpFvPatchScalarField_H

#include "jumpCyclicFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class energyJumpFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class energyJumpFvPatchScalarField
:
    public jumpCyclicFvPatchScalarField
{
    // Private Data

        //- Jump in value from the other patch to this one
        scalarField jump_;


public:

    //- Runtime type information
    TypeName("energyJump");


    // Constructors

        //- Construct from patch and internal field
        energyJumpFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        energyJumpFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given energyJumpFvPatchScalarField onto a
        //  new patch
        energyJumpFvPatchScalarField
        (
            const energyJumpFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        energyJumpFvPatchScalarField
        (
            const energyJumpFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        energyJumpFvPatchScalarField
        (
            const energyJumpFvPatchScalarField&,
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
                new energyJumpFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        // Access

            //- Return the "jump"
            virtual tmp<scalarField> jump() const;


        // Mapping functions

            //- Map the given fvPatchField onto this fvPatchField
            virtual void map(const fvPatchScalarField&, const fieldMapper&);

            //- Reset the fvPatchField to the given fvPatchField
            //  Used for mesh to mesh mapping
            virtual void reset(const fvPatchScalarField&);


        // Evaluation functions

            //- Update the coefficients
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
