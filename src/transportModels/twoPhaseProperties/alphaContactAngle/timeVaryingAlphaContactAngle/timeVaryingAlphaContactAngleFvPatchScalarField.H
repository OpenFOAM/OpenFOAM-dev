/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::timeVaryingAlphaContactAngleFvPatchScalarField

Description
    A time-varying alphaContactAngle scalar boundary condition
    (alphaContactAngleFvPatchScalarField)

SourceFiles
    timeVaryingAlphaContactAngleFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef timeVaryingAlphaContactAngleFvPatchScalarField_H
#define timeVaryingAlphaContactAngleFvPatchScalarField_H

#include "alphaContactAngleFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class timeVaryingAlphaContactAngleFvPatch Declaration
\*---------------------------------------------------------------------------*/

class timeVaryingAlphaContactAngleFvPatchScalarField
:
    public alphaContactAngleFvPatchScalarField
{
    // Private Data

        // Equilibrium contact angle control parameters
        scalar t0_;
        scalar thetaT0_;
        scalar te_;
        scalar thetaTe_;


public:

    //- Runtime type information
    TypeName("timeVaryingAlphaContactAngle");


    // Constructors

        //- Construct from patch and internal field
        timeVaryingAlphaContactAngleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        timeVaryingAlphaContactAngleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  timeVaryingAlphaContactAngleFvPatchScalarField onto a new patch
        timeVaryingAlphaContactAngleFvPatchScalarField
        (
            const timeVaryingAlphaContactAngleFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new timeVaryingAlphaContactAngleFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        timeVaryingAlphaContactAngleFvPatchScalarField
        (
            const timeVaryingAlphaContactAngleFvPatchScalarField&,
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
                new timeVaryingAlphaContactAngleFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        //- Evaluate and return the time-varying equilibrium contact-angle
        virtual tmp<scalarField> theta
        (
            const fvPatchVectorField& Up,
            const fvsPatchVectorField& nHat
        ) const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
