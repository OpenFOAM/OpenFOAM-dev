/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "fanPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* NamedEnum
    <
        fanPressureFvPatchScalarField::fanFlowDirection,
        2
    >::names[] =
    {
        "in",
        "out"
    };
}

const Foam::NamedEnum
<
    Foam::fanPressureFvPatchScalarField::fanFlowDirection,
    2
> Foam::fanPressureFvPatchScalarField::fanFlowDirectionNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fanPressureFvPatchScalarField::fanPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    totalPressureFvPatchScalarField(p, iF),
    fanCurve_(),
    direction_(ffdOut)
{}


Foam::fanPressureFvPatchScalarField::fanPressureFvPatchScalarField
(
    const fanPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    totalPressureFvPatchScalarField(ptf, p, iF, mapper),
    fanCurve_(ptf.fanCurve_),
    direction_(ptf.direction_)
{}


Foam::fanPressureFvPatchScalarField::fanPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    totalPressureFvPatchScalarField(p, iF, dict),
    fanCurve_(dict),
    direction_(fanFlowDirectionNames_.read(dict.lookup("direction")))
{}


Foam::fanPressureFvPatchScalarField::fanPressureFvPatchScalarField
(
    const fanPressureFvPatchScalarField& pfopsf
)
:
    totalPressureFvPatchScalarField(pfopsf),
    fanCurve_(pfopsf.fanCurve_),
    direction_(pfopsf.direction_)
{}


Foam::fanPressureFvPatchScalarField::fanPressureFvPatchScalarField
(
    const fanPressureFvPatchScalarField& pfopsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    totalPressureFvPatchScalarField(pfopsf, iF),
    fanCurve_(pfopsf.fanCurve_),
    direction_(pfopsf.direction_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fanPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Retrieve flux field
    const surfaceScalarField& phi =
        db().lookupObject<surfaceScalarField>(phiName());

    const fvsPatchField<scalar>& phip =
        patch().patchField<surfaceScalarField, scalar>(phi);

    int dir = 2*direction_ - 1;

    // Average volumetric flow rate
    scalar volFlowRate = 0;

    if (phi.dimensions() == dimVelocity*dimArea)
    {
        volFlowRate = dir*gSum(phip);
    }
    else if (phi.dimensions() == dimVelocity*dimArea*dimDensity)
    {
        const scalarField& rhop =
            patch().lookupPatchField<volScalarField, scalar>(rhoName());
        volFlowRate = dir*gSum(phip/rhop);
    }
    else
    {
        FatalErrorIn("fanPressureFvPatchScalarField::updateCoeffs()")
            << "dimensions of phi are not correct"
                << "\n    on patch " << patch().name()
                << " of field " << dimensionedInternalField().name()
                << " in file " << dimensionedInternalField().objectPath() << nl
                << exit(FatalError);
    }

    // Pressure drop for this flow rate
    const scalar pdFan = fanCurve_(max(volFlowRate, 0.0));

    totalPressureFvPatchScalarField::updateCoeffs
    (
        p0() - dir*pdFan,
        patch().lookupPatchField<volVectorField, vector>(UName())
    );
}


void Foam::fanPressureFvPatchScalarField::write(Ostream& os) const
{
    totalPressureFvPatchScalarField::write(os);
    fanCurve_.write(os);
    os.writeKeyword("direction")
        << fanFlowDirectionNames_[direction_] << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fanPressureFvPatchScalarField
    );
};


// ************************************************************************* //
