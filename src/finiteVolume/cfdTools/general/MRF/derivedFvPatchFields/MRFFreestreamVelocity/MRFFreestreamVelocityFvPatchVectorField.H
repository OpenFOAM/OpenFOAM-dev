/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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
    Foam::MRFFreestreamVelocityFvPatchVectorField

Description
    Freestream velocity condition to be used for a wall rotating with the
    moving frame in an MRF (multi-reference frame) or SRF (single reference
    frame) case.

    Given the free stream velocity in the absolute frame at time 0, the
    condition applies the appropriate rotation transformation in time and space
    to determine the current freestream velocity.

Usage
    \table
        Property         | Description                | Required | Default value
        freestreamValue0 | freestream velocity at t=0 | yes      |
        freestreamValue  | freestream velocity at t   | yes      |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type                MRFFreestreamVelocity;
        freestreamValue0    (1 0 0);
        freestreamValue     uniform (1 0 0);
    }
    \endverbatim

See also
    Foam::freestreamFvPatchField
    Foam::MRFVelocityFvPatchVectorField

SourceFiles
    MRFFreestreamVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef MRFFreestreamVelocityFvPatchVectorField_H
#define MRFFreestreamVelocityFvPatchVectorField_H

#include "freestreamVelocityFvPatchVectorField.H"
#include "MRFPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class MRFFreestreamVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class MRFFreestreamVelocityFvPatchVectorField
:
    public freestreamVelocityFvPatchVectorField,
    public MRFPatchField
{
    // Private Data

        //- Velocity of the free stream at time = 0
        vector freestreamValue0_;


public:

    //- Runtime type information
    TypeName("MRFFreestreamVelocity");


    // Constructors

        //- Construct from patch, internal field and dictionary
        MRFFreestreamVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given MRFFreestreamVelocityFvPatchVectorField
        //  onto a new patch
        MRFFreestreamVelocityFvPatchVectorField
        (
            const MRFFreestreamVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        MRFFreestreamVelocityFvPatchVectorField
        (
            const MRFFreestreamVelocityFvPatchVectorField&
        ) = delete;

        //- Copy constructor setting internal field reference
        MRFFreestreamVelocityFvPatchVectorField
        (
            const MRFFreestreamVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new MRFFreestreamVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member Functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
