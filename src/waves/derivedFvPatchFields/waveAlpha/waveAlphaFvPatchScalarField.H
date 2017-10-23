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
    Foam::waveAlphaFvPatchScalarField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides a waveAlpha condition. This sets the phase
    fraction to that specified by a superposition of wave models. All the
    parameters are looked up from the corresponding velocity condition.

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
        inletOutlet true;
    }
    \endverbatim

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
    // Private data

        //- Name of the velocity field
        word UName_;

        //- Is this alpha field that of the liquid under the wave?
        const Switch liquid_;

        //- Act as an inlet/outlet patch? If false, the alpha field will be set
        //  by the wave model even during outflow. This may apply the wave model
        //  more accurately, but it might also be unstable. If true, the alpha
        //  boundary condition will be zero gradient when the flow reverses, as
        //  is usual.
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

        //- Construct as copy
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

        //- Construct as copy setting internal field reference
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


    // Member functions

        // Access

            //- Access the liquid flag
            bool liquid() const
            {
                return liquid_;
            }


        // Evaluation functions

            //- Return the current modelled phase fraction field
            tmp<scalarField> alpha() const;

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
