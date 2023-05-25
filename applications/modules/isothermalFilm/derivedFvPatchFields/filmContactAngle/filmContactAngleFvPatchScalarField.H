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
    Foam::filmContactAngleFvPatchScalarField

Description
    Film contact angle boundary condition.

Usage
    Examples of the contact angle specification:

    Simple constant equilibrium contact angle:
    \verbatim
    <patchName>
    {
        type            filmContactAngle;

        contactAngle
        {
            type                constant;
            theta0              70;
        }
    }
    \endverbatim

    Same as above but using the simpler specification:
    \verbatim
    <patchName>
    {
        type            filmContactAngle;

        theta0          70;
    }
    \endverbatim

    Dynamic contact angle evaluated from the velocity:
    \verbatim
    <patchName>
    {
        type            filmContactAngle;

        contactAngle
        {
            type           dynamic;
            uTheta         1;
            theta0         70;
            thetaAdv       100;
            thetaRec       50;
        }
    }
    \endverbatim


See also
    Foam::contactAngleModel
    Foam::contactAngleModels::constant
    Foam::contactAngleModels::dynamic

SourceFiles
    filmContactAngleFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef filmContactAngleFvPatchScalarField_H
#define filmContactAngleFvPatchScalarField_H

#include "zeroGradientFvPatchFields.H"
#include "contactAngleModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class filmContactAngleFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class filmContactAngleFvPatchScalarField
:
    public zeroGradientFvPatchScalarField
{
    // Private Data

        //- Equilibrium contact angle function
        autoPtr<contactAngleModel> contactAngle_;


public:

    //- Runtime type information
    TypeName("filmContactAngle");


    // Constructors

        //- Construct from patch and internal field
        filmContactAngleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        filmContactAngleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given filmContactAngleFvPatchScalarField
        //  onto a new patch
        filmContactAngleFvPatchScalarField
        (
            const filmContactAngleFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        filmContactAngleFvPatchScalarField
        (
            const filmContactAngleFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        filmContactAngleFvPatchScalarField
        (
            const filmContactAngleFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );


    // Member Functions

        //- Return the cosine of the contact angle
        tmp<scalarField> cosTheta
        (
            const fvPatchVectorField& Up,
            const vectorField& nHat
        ) const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
