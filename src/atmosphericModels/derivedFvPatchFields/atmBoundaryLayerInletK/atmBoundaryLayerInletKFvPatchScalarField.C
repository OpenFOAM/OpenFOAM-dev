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

#include "atmBoundaryLayerInletKFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

atmBoundaryLayerInletKFvPatchScalarField::
atmBoundaryLayerInletKFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(p, iF),
    atmBoundaryLayer()
{}


atmBoundaryLayerInletKFvPatchScalarField::
atmBoundaryLayerInletKFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    inletOutletFvPatchScalarField(p, iF),
    atmBoundaryLayer(patch().Cf(), dict)
{
    phiName_ = dict.lookupOrDefault<word>("phi", "phi");

    refValue() = k(patch().Cf());
    refGrad() = 0;
    valueFraction() = 1;

    if (dict.found("value"))
    {
        scalarField::operator=(scalarField("value", dict, p.size()));
    }
    else
    {
        scalarField::operator=(refValue());
    }
}


atmBoundaryLayerInletKFvPatchScalarField::
atmBoundaryLayerInletKFvPatchScalarField
(
    const atmBoundaryLayerInletKFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    inletOutletFvPatchScalarField(psf, p, iF, mapper),
    atmBoundaryLayer(psf, mapper)
{}


atmBoundaryLayerInletKFvPatchScalarField::
atmBoundaryLayerInletKFvPatchScalarField
(
    const atmBoundaryLayerInletKFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(psf, iF),
    atmBoundaryLayer(psf)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void atmBoundaryLayerInletKFvPatchScalarField::map
(
    const fvPatchScalarField& psf,
    const fvPatchFieldMapper& mapper
)
{
    inletOutletFvPatchScalarField::map(psf, mapper);

    const atmBoundaryLayerInletKFvPatchScalarField& blpsf =
        refCast<const atmBoundaryLayerInletKFvPatchScalarField>(psf);

    atmBoundaryLayer::map(blpsf, mapper);
}


void atmBoundaryLayerInletKFvPatchScalarField::reset
(
    const fvPatchScalarField& psf
)
{
    inletOutletFvPatchScalarField::reset(psf);

    const atmBoundaryLayerInletKFvPatchScalarField& blpsf =
        refCast<const atmBoundaryLayerInletKFvPatchScalarField>(psf);

    atmBoundaryLayer::reset(blpsf);
}


void atmBoundaryLayerInletKFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    atmBoundaryLayer::write(os);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    atmBoundaryLayerInletKFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
