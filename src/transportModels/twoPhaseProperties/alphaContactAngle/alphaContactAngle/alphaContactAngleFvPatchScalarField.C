/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "alphaContactAngleFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(alphaContactAngleFvPatchScalarField, 0);

    template<>
    const char* Foam::NamedEnum
    <
        Foam::alphaContactAngleFvPatchScalarField::limitControls,
        4
    >::names[] =
    {
        "none",
        "gradient",
        "zeroGradient",
        "alpha"
    };
}


const Foam::NamedEnum
<
    Foam::alphaContactAngleFvPatchScalarField::limitControls,
    4
> Foam::alphaContactAngleFvPatchScalarField::limitControlNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::alphaContactAngleFvPatchScalarField::alphaContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    limit_(lcZeroGradient)
{}


Foam::alphaContactAngleFvPatchScalarField::alphaContactAngleFvPatchScalarField
(
    const alphaContactAngleFvPatchScalarField& acpsf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(acpsf, p, iF, mapper),
    limit_(acpsf.limit_)
{}


Foam::alphaContactAngleFvPatchScalarField::alphaContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    limit_(limitControlNames_.read(dict.lookup("limit")))
{
    if (dict.found("gradient"))
    {
        gradient() = scalarField("gradient", dict, p.size());
        fixedGradientFvPatchScalarField::updateCoeffs();
        fixedGradientFvPatchScalarField::evaluate();
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
        gradient() = 0.0;
    }
}


Foam::alphaContactAngleFvPatchScalarField::alphaContactAngleFvPatchScalarField
(
    const alphaContactAngleFvPatchScalarField& acpsf
)
:
    fixedGradientFvPatchScalarField(acpsf),
    limit_(acpsf.limit_)
{}


Foam::alphaContactAngleFvPatchScalarField::alphaContactAngleFvPatchScalarField
(
    const alphaContactAngleFvPatchScalarField& acpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(acpsf, iF),
    limit_(acpsf.limit_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::alphaContactAngleFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{
    if (limit_ == lcGradient)
    {
        gradient() =
        patch().deltaCoeffs()
       *(
           max(min
           (
               *this + gradient()/patch().deltaCoeffs(),
               scalar(1)), scalar(0)
           ) - *this
       );
    }
    else if (limit_ == lcZeroGradient)
    {
        gradient() = 0.0;
    }

    fixedGradientFvPatchScalarField::evaluate();

    if (limit_ == lcAlpha)
    {
        scalarField::operator=(max(min(*this, scalar(1)), scalar(0)));
    }
}


void Foam::alphaContactAngleFvPatchScalarField::write
(
    Ostream& os
) const
{
    fixedGradientFvPatchScalarField::write(os);
    os.writeKeyword("limit")
        << limitControlNames_[limit_] << token::END_STATEMENT << nl;
}


// ************************************************************************* //
