/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

Group
    grpWallBoundaryConditions grpGenericBoundaryConditions

Description
    This boundary condition provides a velocity condition for translational
    motion on walls.

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        U            | translational velocity  | yes|
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            translatingWallVelocity;
        U               (100 0 0);
    }
    \endverbatim


SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    translatingWallVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef translatingWallVelocityFvPatchVectorField_H
#define translatingWallVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

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
    // Private data

        //- Translational velocity
        vector U_;


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

        //- Construct as copy
        translatingWallVelocityFvPatchVectorField
        (
            const translatingWallVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new translatingWallVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
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


    // Member functions

        // Access functions

            //- Return the velocity
            const vector& U() const
            {
                return U_;
            }

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
