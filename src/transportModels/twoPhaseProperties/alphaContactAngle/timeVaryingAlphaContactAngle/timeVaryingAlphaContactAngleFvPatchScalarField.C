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

#include "timeVaryingAlphaContactAngleFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeVaryingAlphaContactAngleFvPatchScalarField::
timeVaryingAlphaContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphaContactAngleFvPatchScalarField(p, iF),
    t0_(0.0),
    thetaT0_(0.0),
    te_(0.0),
    thetaTe_(0.0)
{}


Foam::timeVaryingAlphaContactAngleFvPatchScalarField::
timeVaryingAlphaContactAngleFvPatchScalarField
(
    const timeVaryingAlphaContactAngleFvPatchScalarField& gcpsf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    alphaContactAngleFvPatchScalarField(gcpsf, p, iF, mapper),
    t0_(gcpsf.t0_),
    thetaT0_(gcpsf.thetaT0_),
    te_(gcpsf.te_),
    thetaTe_(gcpsf.te_)
{}


Foam::timeVaryingAlphaContactAngleFvPatchScalarField::
timeVaryingAlphaContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    alphaContactAngleFvPatchScalarField(p, iF, dict),
    t0_(readScalar(dict.lookup("t0"))),
    thetaT0_(readScalar(dict.lookup("thetaT0"))),
    te_(readScalar(dict.lookup("te"))),
    thetaTe_(readScalar(dict.lookup("thetaTe")))
{
    evaluate();
}


Foam::timeVaryingAlphaContactAngleFvPatchScalarField::
timeVaryingAlphaContactAngleFvPatchScalarField
(
    const timeVaryingAlphaContactAngleFvPatchScalarField& gcpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphaContactAngleFvPatchScalarField(gcpsf, iF),
    t0_(gcpsf.t0_),
    thetaT0_(gcpsf.thetaT0_),
    te_(gcpsf.te_),
    thetaTe_(gcpsf.thetaTe_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::timeVaryingAlphaContactAngleFvPatchScalarField::theta
(
    const fvPatchVectorField&,
    const fvsPatchVectorField&
) const
{
    scalar t = patch().boundaryMesh().mesh().time().value();
    scalar theta0 = thetaT0_;

    if (t < t0_)
    {
        theta0 = thetaT0_;
    }
    else if (t > te_)
    {
        theta0 = thetaTe_;
    }
    else
    {
        theta0 = thetaT0_ + (t - t0_)*(thetaTe_ - thetaT0_)/(te_ - t0_);
    }

    return tmp<scalarField>(new scalarField(size(), theta0));
}


void Foam::timeVaryingAlphaContactAngleFvPatchScalarField::write
(
    Ostream& os
) const
{
    alphaContactAngleFvPatchScalarField::write(os);
    os.writeKeyword("t0") << t0_ << token::END_STATEMENT << nl;
    os.writeKeyword("thetaT0") << thetaT0_ << token::END_STATEMENT << nl;
    os.writeKeyword("te") << te_ << token::END_STATEMENT << nl;
    os.writeKeyword("thetaTe") << thetaTe_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        timeVaryingAlphaContactAngleFvPatchScalarField
    );
}

// ************************************************************************* //
