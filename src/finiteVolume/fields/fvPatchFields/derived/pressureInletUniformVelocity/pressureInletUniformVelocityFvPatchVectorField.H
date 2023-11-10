/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::pressureInletUniformVelocityFvPatchVectorField

Description
    This velocity inlet boundary condition is applied to patches where the
    pressure is specified.  The uniform inflow velocity is obtained by
    averaging the flux over the patch, and then applying it in the direction
    normal to the patch faces.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            pressureInletUniformVelocity;
        value           uniform 0;
    }
    \endverbatim

SourceFiles
    pressureInletUniformVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef pressureInletUniformVelocityFvPatchVectorField_H
#define pressureInletUniformVelocityFvPatchVectorField_H

#include "pressureInletVelocityFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
       Class pressureInletUniformVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class pressureInletUniformVelocityFvPatchVectorField
:
    public pressureInletVelocityFvPatchVectorField
{

public:

    //- Runtime type information
    TypeName("pressureInletUniformVelocity");


    // Constructors

        //- Construct from patch, internal field and dictionary
        pressureInletUniformVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  pressureInletUniformVelocityFvPatchVectorField onto a new patch
        pressureInletUniformVelocityFvPatchVectorField
        (
            const pressureInletUniformVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        pressureInletUniformVelocityFvPatchVectorField
        (
            const pressureInletUniformVelocityFvPatchVectorField&
        ) = delete;

        //- Copy constructor setting internal field reference
        pressureInletUniformVelocityFvPatchVectorField
        (
            const pressureInletUniformVelocityFvPatchVectorField&,
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
                new pressureInletUniformVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member Functions

        // Attributes

            //- Return true: this patch field is altered by assignment
            virtual bool assignable() const
            {
                return true;
            }


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();


    // Member Operators

        virtual void operator=(const fvPatchField<vector>& pvf);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
