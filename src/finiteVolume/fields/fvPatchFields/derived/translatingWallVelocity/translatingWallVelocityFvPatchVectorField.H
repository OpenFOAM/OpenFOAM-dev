/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::translatingWallVelocityFvPatchVectorField

Description
    This boundary condition provides a velocity condition for translational
    motion on walls.

Usage
    \table
        Property     | Description             | Required    | Default value
        U            | translational velocity  | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            translatingWallVelocity;
        U               (100 0 0);
    }
    \endverbatim

    The \c U entry is a Function1 of time, see Foam::Function1s.


See also
    Foam::fixedValueFvPatchField
    Foam::Function1s

SourceFiles
    translatingWallVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef translatingWallVelocityFvPatchVectorField_H
#define translatingWallVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class translatingWallVelocityFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class translatingWallVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private Data

        //- Translational velocity
        autoPtr<Function1<vector>> U_;


public:

    //- Runtime type information
    TypeName("translatingWallVelocity");


    // Constructors

        //- Construct from patch and internal field
        translatingWallVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        translatingWallVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given a
        // translatingWallVelocityFvPatchVectorField onto a new patch
        translatingWallVelocityFvPatchVectorField
        (
            const translatingWallVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        translatingWallVelocityFvPatchVectorField
        (
            const translatingWallVelocityFvPatchVectorField&
        ) = delete;

        //- Copy constructor setting internal field reference
        translatingWallVelocityFvPatchVectorField
        (
            const translatingWallVelocityFvPatchVectorField&,
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
                new translatingWallVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member Functions

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
