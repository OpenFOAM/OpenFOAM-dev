/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2023 OpenFOAM Foundation
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

Usage
    \table
        Property  | Description                            | Required? | Default
        phi       | Name of the flux field                 | no        | phi
        liquid    | Is the alpha field that of the liquid? | no        | true
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type        waveAlpha;
        libs        ("libwaves.so");
        liquid      true;
    }
    \endverbatim

See also
    Foam::waveSuperposition

SourceFiles
    waveAlphaFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef waveAlphaFvPatchScalarField_H
#define waveAlphaFvPatchScalarField_H

#include "fixedValueInletOutletFvPatchFields.H"
#include "waveSuperposition.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMeshSubset;

/*---------------------------------------------------------------------------*\
                 Class waveAlphaFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class waveAlphaFvPatchScalarField
:
    public fixedValueInletOutletFvPatchField<scalar>
{
    // Private Data

        //- Is this alpha field that of the liquid under the wave?
        const Switch liquid_;

        //- Mesh subset corresponding to the patch adjacent cells
        mutable autoPtr<fvMeshSubset> faceCellSubset_;

        //- Time index for keeping the subset up to date
        mutable label faceCellSubsetTimeIndex_;


public:

    //- Runtime type information
    TypeName("waveAlpha");


    // Constructors

        //- Construct from patch, internal field and dictionary
        waveAlphaFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedValueTypeFvPatchField onto a new
        //  patch
        waveAlphaFvPatchScalarField
        (
            const waveAlphaFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        waveAlphaFvPatchScalarField
        (
            const waveAlphaFvPatchScalarField&
        ) = delete;

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

            //- Access the face-cell subset
            const fvMeshSubset& faceCellSubset() const;


        // Evaluation functions

            //- Return the current modelled phase fraction field on the patch
            //  faces at the given time
            tmp<scalarField> alpha(const scalar t) const;

            //- Return the current modelled phase fraction field in the
            //  neighbour cells at the given time
            tmp<scalarField> alphan(const scalar t) const;

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
