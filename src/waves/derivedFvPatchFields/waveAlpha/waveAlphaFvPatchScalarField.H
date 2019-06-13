/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2019 OpenFOAM Foundation
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
    Foam::waveAlphaFvPatchScalarField

Description
    This boundary condition provides a waveAlpha condition. This sets the phase
    fraction to that specified by a superposition of wave models. All the
    wave modelling parameters are obtained from a centrally registered
    waveSuperposition class.

    Flow reversal will occur in the event that the amplitude of the velocity
    oscillation is greater than the mean flow. This triggers special handling,
    the form of which depends on the inletOutlet flag and whether a wave
    pressure condition is being used.

    If a wave pressure condition is not being used, the inletOutlet switches
    between a fixedValue and an inletOutlet condition, with the value given by
    the wave model. If fixedValue, the result may be more accurate, but it
    might also be unstable.

    If a wave pressure condition is being used, then the normal phase fraction
    condition becomes fixedGradient on outlet faces. This gradient is
    calculated numerically by evaluating the wave model on both the patch face
    and the adjacent cell.

Usage
    \table
        Property  | Description             | Req'd? | Default
        U         | name of the velocity field | no  | U
        liquid    | is the alpha field that of the liquid | no | true
        inletOutlet | does the condition behave like inletOutlet | no | true
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type        waveAlpha;
        U           U;
        liquid      true;
        inletOutlet true;
    }
    \endverbatim

See also
    Foam::waveSuperposition

SourceFiles
    waveAlphaFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef waveAlphaFvPatchScalarField_H
#define waveAlphaFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class waveAlphaFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class waveAlphaFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private Data

        //- Name of the velocity field
        const word UName_;

        //- Is this alpha field that of the liquid under the wave?
        const Switch liquid_;

        //- Act as an inlet/outlet patch?
        const Switch inletOutlet_;


public:

    //- Runtime type information
    TypeName("waveAlpha");


    // Constructors

        //- Construct from patch and internal field
        waveAlphaFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        waveAlphaFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given mixedTypeFvPatchField onto a new patch
        waveAlphaFvPatchScalarField
        (
            const waveAlphaFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        waveAlphaFvPatchScalarField
        (
            const waveAlphaFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new waveAlphaFvPatchScalarField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        waveAlphaFvPatchScalarField
        (
            const waveAlphaFvPatchScalarField&,
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
                new waveAlphaFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        // Access

            //- Access the liquid flag
            bool liquid() const
            {
                return liquid_;
            }


        // Evaluation functions

            //- Return the current modelled phase fraction field
            tmp<scalarField> alpha() const;

            //- Return the current modelled phase fraction field in the
            //  neighbour cell
            tmp<scalarField> alphan() const;

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
