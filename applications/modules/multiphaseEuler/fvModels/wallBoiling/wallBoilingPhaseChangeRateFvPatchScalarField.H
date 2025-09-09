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
    Foam::wallBoilingPhaseChangeRateFvPatchScalarField

Description
    This boundary condition is used for the phase change rate field of the wall
    boiling fvModel. It contains the phase change rate, and also all the other
    state fields used by the wall boiling model. The phase change rate field
    and its boundary conditions are constructed automatically so the user need
    never manually specify this condition.

See also
    Foam::fv::wallBoiling

SourceFiles
    wallBoilingPhaseChangeRateFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef wallBoilingPhaseChangeRateFvPatchScalarField_H
#define wallBoilingPhaseChangeRateFvPatchScalarField_H

#include "wallBoiling.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class wallBoilingPhaseChangeRateFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class wallBoilingPhaseChangeRateFvPatchScalarField
:
    public calculatedFvPatchScalarField
{
    // Private Static Data

        //- Name-to-index map for the property names
        static const hashedWordList propertyNames_;

        //- Property pointers
        static const List
        <
            const scalarField wallBoilingPhaseChangeRateFvPatchScalarField::*
        > propertyPtrs_;

        //- Property dimensions
        static const List<const dimensionSet*> propertyDimensions_;


    // Private Data

        //- Wall liquid fraction
        scalarField wetFraction_;

        //- Departure diameter
        scalarField dDeparture_;

        //- Departure frequency
        scalarField fDeparture_;

        //- Nucleation site density
        scalarField nucleationSiteDensity_;

        //- Quenching surface heat flux
        scalarField qQuenching_;

        //- Evaporative surface heat flux
        scalarField qEvaporative_;

        //- Liquid turbulent thermal diffusivity
        scalarField alphatLiquid_;

        //- Vapour turbulent thermal diffusivity
        scalarField alphatVapour_;


public:

    //- Allow the wall boiling model to access and modify the state fields
    friend class fv::wallBoiling;


    //- Runtime type information
    TypeName("wallBoilingPhaseChangeRate");


    // Constructors

        //- Construct from patch and internal field
        wallBoilingPhaseChangeRateFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        wallBoilingPhaseChangeRateFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  wallBoilingPhaseChangeRateFvPatchScalarField
        //  onto a new patch
        wallBoilingPhaseChangeRateFvPatchScalarField
        (
            const wallBoilingPhaseChangeRateFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        wallBoilingPhaseChangeRateFvPatchScalarField
        (
            const wallBoilingPhaseChangeRateFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        wallBoilingPhaseChangeRateFvPatchScalarField
        (
            const wallBoilingPhaseChangeRateFvPatchScalarField&,
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
                new wallBoilingPhaseChangeRateFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        //- Access one of the property fields by name
        const scalarField& property(const word& name) const;

        //- Access one of the property fields' dimensions by name
        static const dimensionSet& propertyDimensions(const word& name);

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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
