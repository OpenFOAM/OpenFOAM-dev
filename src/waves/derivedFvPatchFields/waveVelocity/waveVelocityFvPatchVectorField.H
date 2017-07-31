/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
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
    Foam::waveVelocityFvPatchVectorField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides a waveVelocity condition. This sets the
    velocity to that specified by a superposition of wave models. The
    corresponding phase fraction condition looks this condition up and re-uses
    the wave modelling.

Usage
    \table
        Property    | Description                            | Req'd? | Default
        origin      | origin of the wave coordinate system   | yes    |
        direction   | direction of the mean flow             | yes    |
        speed       | speed of the mean flow                 | yes    |
        waves       | list of wave models to superimpose     | yes    |
        scale       | scale factor along the mean flow direction | no | None
        crossScale  | scale factor across the mean flow direction | no | None
        phi         | Name of the flux field                 | no     | phi
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type        waveVelocity;
        origin      (0 25 0);
        direction   (1 0 0);
        speed       2;
        waves
        (
            Airy
            {
                length      40;
                amplitude   0.5;
                phase       0;
                angle       0;
            }
            Airy
            {
                length      20;
                amplitude   0.25;
                phase       1.5708;
                angle       0;
            }
        );
        scale       table ((100 1) (200 0));
        crossScale  constant 1;
    }
    \endverbatim

SourceFiles
    waveVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef waveVelocityFvPatchVectorField_H
#define waveVelocityFvPatchVectorField_H

#include "directionMixedFvPatchFields.H"
#include "waveSuperposition.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class waveVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class waveVelocityFvPatchVectorField
:
    public directionMixedFvPatchVectorField
{
    // Private data

        //- Name of the flux field
        const word phiName_;

        //- Wave superposition
        const waveSuperposition waves_;


public:

    //- Runtime type information
    TypeName("waveVelocity");


    // Constructors

        //- Construct from patch and internal field
        waveVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        waveVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given mixedTypeFvPatchField onto a new patch
        waveVelocityFvPatchVectorField
        (
            const waveVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        waveVelocityFvPatchVectorField
        (
            const waveVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new waveVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        waveVelocityFvPatchVectorField
        (
            const waveVelocityFvPatchVectorField&,
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
                new waveVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Access the wave models
            const waveSuperposition& waves() const
            {
                return waves_;
            }


        // Evaluation functions

            //- Return the current modelled velocity field
            tmp<vectorField> U() const;

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
