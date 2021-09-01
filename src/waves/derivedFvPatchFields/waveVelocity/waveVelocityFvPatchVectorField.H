/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2021 OpenFOAM Foundation
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

Description
    This boundary condition provides a waveVelocity condition. This sets the
    velocity to that specified by a superposition of wave models. All the
    wave modelling parameters are obtained from a centrally registered
    waveSuperposition class.

Usage
    \table
        Property  | Description                            | Req'd? | Default
        phi       | Name of the flux field                 | no     | phi
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type        waveVelocity;
        phi         phi;
    }
    \endverbatim

See also
    Foam::waveSuperposition

SourceFiles
    waveVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef waveVelocityFvPatchVectorField_H
#define waveVelocityFvPatchVectorField_H

#include "fixedValueInletOutletFvPatchFields.H"
#include "waveSuperposition.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMeshSubset;

/*---------------------------------------------------------------------------*\
               Class waveVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class waveVelocityFvPatchVectorField
:
    public fixedValueInletOutletFvPatchField<vector>
{
    // Private Data

        //- Mesh subset corresponding to the patch adjacent cells
        mutable autoPtr<fvMeshSubset> faceCellSubset_;

        //- Time index for keeping the subset up to date
        mutable label faceCellSubsetTimeIndex_;


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

        //- Disallow copy without setting internal field reference
        waveVelocityFvPatchVectorField
        (
            const waveVelocityFvPatchVectorField&
        ) = delete;

        //- Copy constructor setting internal field reference
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


    // Member Functions

        // Access

            //- Access the face-cell subset
            const fvMeshSubset& faceCellSubset() const;


        // Evaluation functions

            //- Return the current modelled velocity field on the patch faces
            //  at the given time
            tmp<vectorField> U(const scalar t) const;

            //- Return the current modelled velocity field in the neighbour
            //  cells at the given time
            tmp<vectorField> Un(const scalar t) const;

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
