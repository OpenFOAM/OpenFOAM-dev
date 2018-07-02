/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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
    Foam::fixedMultiPhaseHeatFluxFvPatchScalarField

Description
    Calculates a wall temperature that produces the specified overall wall heat
    flux across all the phases in an Eulerian multi-phase simulation.

    Intended to be used with copiedFixedValue to ensure that phase wall
    temperature are consistent:
        - Set 'fixedMultiPhaseHeatFlux' boundary for one of the phases
        - Use 'copiedFixedValue' for all the other phases.

See also
    Foam::fixedValueFvPatchField

SourceFiles
    fixedMultiPhaseHeatFluxFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedMultiPhaseHeatFluxFvPatchScalarField_H
#define fixedMultiPhaseHeatFluxFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class fixedMultiPhaseHeatFluxFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class fixedMultiPhaseHeatFluxFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Heat power [W] or flux [W/m2]
        scalarField q_;

        //- Relaxation factor
        scalar relax_;

        //- Minimum temperature limit [K]
        scalar Tmin_;


public:

    //- Runtime type information
    TypeName("fixedMultiPhaseHeatFlux");


    // Constructors

        //- Construct from patch and internal field
        fixedMultiPhaseHeatFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedMultiPhaseHeatFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  fixedMultiPhaseHeatFluxFvPatchScalarField
        //  onto a new patch
        fixedMultiPhaseHeatFluxFvPatchScalarField
        (
            const fixedMultiPhaseHeatFluxFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixedMultiPhaseHeatFluxFvPatchScalarField
        (
            const fixedMultiPhaseHeatFluxFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        fixedMultiPhaseHeatFluxFvPatchScalarField
        (
            const fixedMultiPhaseHeatFluxFvPatchScalarField&,
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
