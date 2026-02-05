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
    Foam::wallCondensationPhaseChangeRateFvPatchScalarField

Description
    This boundary condition is used for the phase change rate field of the wall
    condensation fvModel. It contains the phase change rate, and also all the
    other state fields used by the wall condensation model. The phase change
    rate field and its boundary conditions are constructed automatically so
    the user need never manually specify this condition.

See also
    Foam::fv::wallCondensation

SourceFiles
    wallCondensationPhaseChangeRateFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef wallCondensationPhaseChangeRateFvPatchScalarField_H
#define wallCondensationPhaseChangeRateFvPatchScalarField_H

#include "wallCondensation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class wallCondensationPhaseChangeRateFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class wallCondensationPhaseChangeRateFvPatchScalarField
:
    public calculatedFvPatchScalarField
{
    // Private Data

        //- Mask which indicates which faces are condensing
        scalarField condensing_;

        //- Liquid turbulent thermal diffusivity
        scalarField alphatLiquid_;

        //- Vapour turbulent thermal diffusivity
        scalarField alphatVapour_;


public:

    //- Allow the wall condensation model to access and modify the state fields
    friend class fv::wallCondensation;


    //- Runtime type information
    TypeName("wallCondensationPhaseChangeRate");


    // Constructors

        //- Construct from patch and internal field
        wallCondensationPhaseChangeRateFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        wallCondensationPhaseChangeRateFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  wallCondensationPhaseChangeRateFvPatchScalarField
        //  onto a new patch
        wallCondensationPhaseChangeRateFvPatchScalarField
        (
            const wallCondensationPhaseChangeRateFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        wallCondensationPhaseChangeRateFvPatchScalarField
        (
            const wallCondensationPhaseChangeRateFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        wallCondensationPhaseChangeRateFvPatchScalarField
        (
            const wallCondensationPhaseChangeRateFvPatchScalarField&,
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
                new wallCondensationPhaseChangeRateFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        //- Access the liquid turbulent thermal diffusivity
        const scalarField& alphatLiquid() const;

        //- Access the vapour turbulent thermal diffusivity
        const scalarField& alphatVapour() const;


        // Mapping functions

            //- Map the given fvPatchField onto this fvPatchField
            virtual void map(const fvPatchScalarField&, const fieldMapper&);

            //- Reset the fvPatchField to the given fvPatchField
            //  Used for mesh to mesh mapping
            virtual void reset(const fvPatchScalarField&);


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            virtual void write(Ostream&) const;


    // Member Operators

        //- Inherit assignment
        using calculatedFvPatchScalarField::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
