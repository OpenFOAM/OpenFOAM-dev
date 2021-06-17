/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "dynamicAlphaContactAngleFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicAlphaContactAngleFvPatchScalarField::
dynamicAlphaContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphaContactAngleFvPatchScalarField(p, iF),
    theta0_(0.0),
    uTheta_(0.0),
    thetaA_(0.0),
    thetaR_(0.0)
{}


Foam::dynamicAlphaContactAngleFvPatchScalarField::
dynamicAlphaContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    alphaContactAngleFvPatchScalarField(p, iF, dict),
    theta0_(dict.lookup<scalar>("theta0")),
    uTheta_(dict.lookup<scalar>("uTheta")),
    thetaA_(dict.lookupBackwardsCompatible<scalar>({"thetaA", "thetaRec"})),
    thetaR_(dict.lookupBackwardsCompatible<scalar>({"thetaR", "thetaAdv"}))
{
    evaluate();
}


Foam::dynamicAlphaContactAngleFvPatchScalarField::
dynamicAlphaContactAngleFvPatchScalarField
(
    const dynamicAlphaContactAngleFvPatchScalarField& gcpsf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    alphaContactAngleFvPatchScalarField(gcpsf, p, iF, mapper),
    theta0_(gcpsf.theta0_),
    uTheta_(gcpsf.uTheta_),
    thetaA_(gcpsf.thetaA_),
    thetaR_(gcpsf.thetaR_)
{}


Foam::dynamicAlphaContactAngleFvPatchScalarField::
dynamicAlphaContactAngleFvPatchScalarField
(
    const dynamicAlphaContactAngleFvPatchScalarField& gcpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphaContactAngleFvPatchScalarField(gcpsf, iF),
    theta0_(gcpsf.theta0_),
    uTheta_(gcpsf.uTheta_),
    thetaA_(gcpsf.thetaA_),
    thetaR_(gcpsf.thetaR_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::dynamicAlphaContactAngleFvPatchScalarField::theta
(
    const fvPatchVectorField& Up,
    const fvsPatchVectorField& nHat
) const
{
    if (uTheta_ < small)
    {
        return tmp<scalarField>(new scalarField(size(), theta0_));
    }

    const vectorField nf(patch().nf());

    // Calculated the component of the velocity parallel to the wall
    vectorField Uwall(Up.patchInternalField() - Up);
    Uwall -= (nf & Uwall)*nf;

    // Find the direction of the interface parallel to the wall
    vectorField nWall(nHat - (nf & nHat)*nf);

    // Normalise nWall
    nWall /= (mag(nWall) + small);

    // Calculate Uwall resolved normal to the interface parallel to
    // the interface
    scalarField uwall(nWall & Uwall);

    return theta0_ + (thetaA_ - thetaR_)*tanh(uwall/uTheta_);
}


void Foam::dynamicAlphaContactAngleFvPatchScalarField::write(Ostream& os) const
{
    alphaContactAngleFvPatchScalarField::write(os);
    writeEntry(os, "theta0", theta0_);
    writeEntry(os, "uTheta", uTheta_);
    writeEntry(os, "thetaA", thetaA_);
    writeEntry(os, "thetaR", thetaR_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        dynamicAlphaContactAngleFvPatchScalarField
    );
}


// ************************************************************************* //
