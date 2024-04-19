/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    totalPressureFvPatchScalarField(p, iF, dict),
    fanCurve_
    (
        Function1<scalar>::New
        (
            "fanCurve",
            dimVolumetricFlux,
            iF.dimensions(),
            dict
        )
    ),
    direction_(fanFlowDirectionNames_.read(dict.lookup("direction")))
{}


Foam::fanPressureFvPatchScalarField::fanPressureFvPatchScalarField
(
    const fanPressureFvPatchScalarField& pfopsf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    totalPressureFvPatchScalarField(pfopsf, p, iF, mapper),
    fanCurve_(pfopsf.fanCurve_, false),
    direction_(pfopsf.direction_)
{}


Foam::fanPressureFvPatchScalarField::fanPressureFvPatchScalarField
(
    const fanPressureFvPatchScalarField& pfopsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    totalPressureFvPatchScalarField(pfopsf, iF),
    fanCurve_(pfopsf.fanCurve_, false),
    direction_(pfopsf.direction_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fanPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    const fvPatchField<vector>& Up =
        patch().lookupPatchField<volVectorField, vector>(UName_);

    const int sign = direction_ ? +1 : -1;

    // Get the volumetric flow rate
    scalar volFlowRate = 0;
    if (phip.internalField().dimensions() == dimVolumetricFlux)
    {
        volFlowRate = sign*gSum(phip);
    }
    else if
    (
        phip.internalField().dimensions() == dimMassFlux
    )
    {
        const scalarField& rhop =
            patch().lookupPatchField<volScalarField, scalar>(rhoName_);

        volFlowRate = sign*gSum(phip/rhop);
    }
    else
    {
        FatalErrorInFunction
            << "dimensions of phi are not correct"
            << "\n    on patch " << patch().name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath() << nl
            << exit(FatalError);
    }

    // Pressure drop for this flow rate
    const scalar dp0 = fanCurve_->value(max(volFlowRate, scalar(0)));

    dynamicPressureFvPatchScalarField::updateCoeffs
    (
        p0_ - sign*dp0,
        -0.5*neg(phip)*magSqr(Up)
    );
}


void Foam::fanPressureFvPatchScalarField::write(Ostream& os) const
{
    totalPressureFvPatchScalarField::write(os);
    writeEntry(os, fanCurve_());
    writeEntry(os, "direction", fanFlowDirectionNames_[direction_]);
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
