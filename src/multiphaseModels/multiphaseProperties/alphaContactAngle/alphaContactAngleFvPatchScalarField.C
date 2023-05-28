/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::alphaContactAngleFvPatchScalarField::alphaContactAngleFvPatchScalarField
(
    const alphaContactAngleFvPatchScalarField& gcpsf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    zeroGradientFvPatchScalarField(gcpsf, p, iF, mapper),
    thetaProps_(gcpsf.thetaProps_)
{}


Foam::alphaContactAngleFvPatchScalarField::alphaContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    zeroGradientFvPatchScalarField(p, iF, dict),
    thetaProps_()
{
    forAllConstIter(dictionary, dict.subDict("contactAngleProperties"), iter)
    {
        thetaProps_.insert
        (
            iter().keyword(),
            contactAngleProperties(iter().dict())
        );
    }
}


Foam::alphaContactAngleFvPatchScalarField::alphaContactAngleFvPatchScalarField
(
    const alphaContactAngleFvPatchScalarField& gcpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    zeroGradientFvPatchScalarField(gcpsf, iF),
    thetaProps_(gcpsf.thetaProps_)
{}


Foam::alphaContactAngleFvPatchScalarField::alphaContactAngleFvPatchScalarField::
contactAngleProperties::contactAngleProperties()
:
    theta0_(NaN),
    dynamic_(false),
    uTheta_(NaN),
    thetaA_(NaN),
    thetaR_(NaN)
{}


Foam::alphaContactAngleFvPatchScalarField::alphaContactAngleFvPatchScalarField::
contactAngleProperties::contactAngleProperties(const scalar theta0)
:
    theta0_(theta0),
    dynamic_(false),
    uTheta_(NaN),
    thetaA_(NaN),
    thetaR_(NaN)
{}


Foam::alphaContactAngleFvPatchScalarField::alphaContactAngleFvPatchScalarField::
contactAngleProperties::contactAngleProperties
(
    const scalar theta0,
    const scalar uTheta,
    const scalar thetaA,
    const scalar thetaR
)
:
    theta0_(theta0),
    dynamic_(true),
    uTheta_(uTheta),
    thetaA_(thetaA),
    thetaR_(thetaR)
{}

Foam::alphaContactAngleFvPatchScalarField::alphaContactAngleFvPatchScalarField::
contactAngleProperties::contactAngleProperties(const dictionary& dict)
:
    theta0_(dict.lookup<scalar>("theta0")),
    dynamic_(dict.found("uTheta")),
    uTheta_(dynamic_ ? dict.lookup<scalar>("uTheta") : NaN),
    thetaA_(dynamic_ ? dict.lookup<scalar>("thetaA") : NaN),
    thetaR_(dynamic_ ? dict.lookup<scalar>("thetaR") : NaN)
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::alphaContactAngleFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);

    writeKeyword(os, "contactAngleProperties")
        << nl << indent << token::BEGIN_BLOCK << nl << incrIndent;
    forAllConstIter(HashTable<contactAngleProperties>, thetaProps_, iter)
    {
        writeKeyword(os, iter.key())
            << nl << indent << token::BEGIN_BLOCK << nl << incrIndent;
        iter().write(os);
        os << decrIndent << indent << token::END_BLOCK << endl;
    }
    os << decrIndent << indent << token::END_BLOCK << endl;

    writeEntry(os, "value", *this);
}


Foam::alphaContactAngleFvPatchScalarField::alphaContactAngleFvPatchScalarField::
contactAngleProperties
Foam::alphaContactAngleFvPatchScalarField::alphaContactAngleFvPatchScalarField::
contactAngleProperties::reversed()
const
{
    return
        dynamic()
      ? contactAngleProperties
        (
            180 - theta0_,
            uTheta_,
            180 - thetaA_,
            180 - thetaR_
        )
      : contactAngleProperties
        (
            180 - theta0_
        );
}


void
Foam::alphaContactAngleFvPatchScalarField::alphaContactAngleFvPatchScalarField::
contactAngleProperties::write(Ostream& os) const
{
    writeEntry(os, "theta0", theta0_);
    if (dynamic())
    {
        writeEntry(os, "uTheta", uTheta_);
        writeEntry(os, "thetaA", thetaA_);
        writeEntry(os, "thetaR", thetaR_);
    }
}


bool
Foam::alphaContactAngleFvPatchScalarField::alphaContactAngleFvPatchScalarField::
contactAngleProperties::operator==
(
    const contactAngleProperties& thetaProps
) const
{
    if (dynamic() != thetaProps.dynamic()) return false;

    static const scalar thetaTol = 180*rootSmall;
    static const scalar uThetaTol = rootSmall;

    return
        dynamic()
      ? mag(theta0() - thetaProps.theta0()) < thetaTol
     && mag(uTheta() - thetaProps.uTheta()) < uThetaTol
     && mag(thetaA() - thetaProps.thetaA()) < thetaTol
     && mag(thetaR() - thetaProps.thetaR()) < thetaTol
      : mag(theta0() - thetaProps.theta0()) < thetaTol;
}


bool
Foam::alphaContactAngleFvPatchScalarField::alphaContactAngleFvPatchScalarField::
contactAngleProperties::operator!=
(
    const contactAngleProperties& thetaProps
) const
{
    return !(*this == thetaProps);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        alphaContactAngleFvPatchScalarField
    );
}

// ************************************************************************* //
