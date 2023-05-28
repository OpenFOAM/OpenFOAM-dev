/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::filmSurfaceVelocityFvPatchVectorField

Description
    Film surface velocity boundary condition

    Evaluates the surface velocity from the shear imposed by the neighbouring
    fluid velocity using either a simple drag model based on the difference
    between the fluid and film velocities multiplied by the coefficient \c Cs or
    if \c Cs is not specified or set to 0 the fluid viscous shear stress.

    The simple model might be used in preference to the fluid viscous shear
    stress model in order to provide some means to include the drag enhancing
    effect of surface ripples, rivulets etc. in the film surface.

Usage
    \table
        Property     | Description             | Required    | Default value
        Cs           | Fluid-film drag coefficient | no | 0
    \endtable

    Example of the boundary condition specification using the simple drag model:
    \verbatim
    <patchName>
    {
        type            filmSurfaceVelocity;
        Cs              0.005;
        value           $internalField;
    }
    \endverbatim

    Example of the boundary condition specification using the fluid stress:
    \verbatim
    <patchName>
    {
        type            filmSurfaceVelocity;
        value           $internalField;
    }
    \endverbatim

See also
    Foam::mixedFvPatchField

SourceFiles
    filmSurfaceVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef filmSurfaceVelocityFvPatchVectorField_H
#define filmSurfaceVelocityFvPatchVectorField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
          Class filmSurfaceVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class filmSurfaceVelocityFvPatchVectorField
:
    public mixedFvPatchVectorField
{
    // Private Data

        //- Fluid-film drag-coefficient
        scalar Cs_;


public:

   //- Runtime type information
   TypeName("filmSurfaceVelocity");


   // Constructors

        //- Construct from patch, internal field and dictionary
        filmSurfaceVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  filmSurfaceVelocityFvPatchVectorField
        //  onto a new patch
        filmSurfaceVelocityFvPatchVectorField
        (
            const filmSurfaceVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        filmSurfaceVelocityFvPatchVectorField
        (
            const filmSurfaceVelocityFvPatchVectorField&
        ) = delete;

        //- Copy constructor setting internal field reference
        filmSurfaceVelocityFvPatchVectorField
        (
            const filmSurfaceVelocityFvPatchVectorField&,
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
                new filmSurfaceVelocityFvPatchVectorField
                (
                    *this,
                    iF
                )
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
