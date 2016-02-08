/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 OpenFOAM Foundation
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
    Foam::temperatureDependentAlphaContactAngleFvPatchScalarField

Description
    Temperature-dependent constant alphaContactAngle scalar boundary condition.

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        TName        | Temperature field name  | no          | T
        theta0       | Contact angle data      | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            temperatureDependentAlphaContactAngle;
        theta0          constant 60;
    }
    \endverbatim

SeeAlso
    Foam::alphaContactAngleFvPatchScalarField
    Foam::constantAlphaContactAngleFvPatchScalarField

SourceFiles
    temperatureDependentAlphaContactAngleFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef temperatureDependentAlphaContactAngleFvPatchScalarField_H
#define temperatureDependentAlphaContactAngleFvPatchScalarField_H

#include "alphaContactAngleFvPatchScalarField.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
       Class temperatureDependentAlphaContactAngleFvPatch Declaration
\*---------------------------------------------------------------------------*/

class temperatureDependentAlphaContactAngleFvPatchScalarField
:
    public alphaContactAngleFvPatchScalarField
{
    // Private data

        //- Name of temperature field, default = "T"
        word TName_;

        //- Equilibrium contact angle table
        autoPtr<Function1<scalar>> theta0_;


public:

    //- Runtime type information
    TypeName("temperatureDependentAlphaContactAngle");


    // Constructors

        //- Construct from patch and internal field
        temperatureDependentAlphaContactAngleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        temperatureDependentAlphaContactAngleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  temperatureDependentAlphaContactAngleFvPatchScalarField
        //  onto a new patch
        temperatureDependentAlphaContactAngleFvPatchScalarField
        (
            const temperatureDependentAlphaContactAngleFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        temperatureDependentAlphaContactAngleFvPatchScalarField
        (
            const temperatureDependentAlphaContactAngleFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new temperatureDependentAlphaContactAngleFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        temperatureDependentAlphaContactAngleFvPatchScalarField
        (
            const temperatureDependentAlphaContactAngleFvPatchScalarField&,
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
                new temperatureDependentAlphaContactAngleFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        //- Return the equilibrium contact-angle
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
