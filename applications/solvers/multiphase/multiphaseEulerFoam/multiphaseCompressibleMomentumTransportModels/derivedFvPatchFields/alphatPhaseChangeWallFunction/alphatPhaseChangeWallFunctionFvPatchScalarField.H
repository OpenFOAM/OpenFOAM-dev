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
    Foam::compressible::alphatPhaseChangeWallFunctionFvPatchScalarField

Description
    Abstract base-class for all alphatWallFunctions supporting phase-change.

See also
    Foam::alphatPhaseJayatillekeWallFunctionFvPatchScalarField

SourceFiles
    alphatPhaseChangeWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef alphatPhaseChangeWallFunctionFvPatchScalarField_H
#define alphatPhaseChangeWallFunctionFvPatchScalarField_H

#include "alphatPhaseJayatillekeWallFunctionFvPatchScalarField.H"
#include "phasePairKey.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
       Class alphatPhaseChangeWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class alphatPhaseChangeWallFunctionFvPatchScalarField
:
    public alphatPhaseJayatillekeWallFunctionFvPatchScalarField
{
protected:

    // Protected data

        //- Name of the other phase
        const word otherPhaseName_;

        //- Relaxation factor
        const scalar relax_;

        //- Rate of phase-change
        scalarField dmdtf_;


public:

    //- Runtime type information
    TypeName("compressible::alphatPhaseChangeWallFunction");


    // Constructors

        //- Construct from patch and internal field
        alphatPhaseChangeWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        alphatPhaseChangeWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  alphatPhaseChangeWallFunctionFvPatchScalarField
        //  onto a new patch
        alphatPhaseChangeWallFunctionFvPatchScalarField
        (
            const alphatPhaseChangeWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        alphatPhaseChangeWallFunctionFvPatchScalarField
        (
            const alphatPhaseChangeWallFunctionFvPatchScalarField&
        );

        //- Copy constructor setting internal field reference
        alphatPhaseChangeWallFunctionFvPatchScalarField
        (
            const alphatPhaseChangeWallFunctionFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );


    // Member Functions

        //- Is there phase change mass transfer for this phasePair
        bool activePhasePair(const phasePairKey&) const;

        //- Return the rate of phase-change
        const scalarField& dmdtf() const;

        //- Return the rate of phase-change for specific phase pair
        const scalarField& dmdtf(const phasePairKey&) const;


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            //  Used to update fields following mesh topology change
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            //  Used to reconstruct fields
            virtual void rmap(const fvPatchScalarField&, const labelList&);


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs() = 0;


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
