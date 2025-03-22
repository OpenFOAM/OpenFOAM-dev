/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024-2025 OpenFOAM Foundation
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

#include "lumpedMassTemperatureFvPatchScalarField.H"
#include "fieldMapper.H"
#include "thermophysicalTransportModel.H"
#include "ZeroConstant.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::lumpedMassTemperatureFvPatchScalarField::closed() const
{
    return mag(gSum(patch().Sf()))/gSum(patch().magSf()) < rootSmall;
}


Foam::scalar Foam::lumpedMassTemperatureFvPatchScalarField::V() const
{
    return
       -gSum(patch().Sf() & patch().Cf())
       /patch().boundaryMesh().mesh().nSolutionD();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lumpedMassTemperatureFvPatchScalarField::
lumpedMassTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    rho_(dict.lookup<scalar>("rho")),
    Cv_(dict.lookup<scalar>("Cv")),
    T_
    (
        IOobject
        (
            "T_" + patch().name(),
            db().time().name(),
            db()
        ),
        dimensionedScalar(dimTemperature, dict.lookup<scalar>("T"))
    ),
    Q_
    (
        dict.found("Q")
      ? Function1<scalar>::New
        (
            "Q",
            db().time().userUnits(),
            dimPower,
            dict
        )
      : autoPtr<Function1<scalar>>(new Function1s::ZeroConstant<scalar>("Q"))
    ),
    V_(NaN)
{
    if (!dict.readIfPresent("volume", V_))
    {
        if (closed())
        {
            V_ = V();
            Info<< "Volume for the thermal mass, enclosed by patch '"
                << patch().name() << "', = " << V_;
        }
        else
        {
            FatalErrorInFunction
                << "Patch '" << patch().name()
                << "' corresponding to a thermal mass is not closed." << nl
                << "Please specify the volume with the optional 'volume' entry."
                << exit(FatalError);
        }
    }

    fvPatchScalarField::operator=(T_.value());
}


Foam::lumpedMassTemperatureFvPatchScalarField::
lumpedMassTemperatureFvPatchScalarField
(
    const lumpedMassTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    rho_(ptf.rho_),
    Cv_(ptf.Cv_),
    T_(ptf.T_),
    Q_(ptf.Q_),
    V_(ptf.V_)
{}


Foam::lumpedMassTemperatureFvPatchScalarField::
lumpedMassTemperatureFvPatchScalarField
(
    const lumpedMassTemperatureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    rho_(ptf.rho_),
    Cv_(ptf.Cv_),
    T_(ptf.T_),
    Q_(ptf.Q_),
    V_(ptf.V_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::lumpedMassTemperatureFvPatchScalarField::map
(
    const fvPatchScalarField& ptf,
    const fieldMapper& mapper
)
{
    fixedValueFvPatchScalarField::map(ptf, mapper);
}


void Foam::lumpedMassTemperatureFvPatchScalarField::reset
(
    const fvPatchScalarField& ptf
)
{
    fixedValueFvPatchScalarField::reset(ptf);
}


void Foam::lumpedMassTemperatureFvPatchScalarField::updateCoeffs()
{
    if
    (
        updated()
    )
    {
        return;
    }

    const thermophysicalTransportModel& ttm =
        db().lookupType<thermophysicalTransportModel>
        (
            internalField().group()
        );

    const scalarField Hf
    (
        ttm.kappaEff(patch().index())*patch().magSf()*patch().deltaCoeffs()
    );

    const scalar Hs = rho_*Cv_*V_/db().time().deltaTValue();

    T_.value() =
    (
        Q_->value(db().time().value())
      + gSum(Hf*patchInternalField())
      + Hs*T_.oldTime().value()
    )/(Hs + gSum(Hf));

    operator==(T_.value());

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::lumpedMassTemperatureFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    writeEntry(os, "rho", rho_);
    writeEntry(os, "Cv", Cv_);
    writeEntry(os, "T", T_.value());
    writeEntry(os, Q_());
    writeEntry(os, "volume", V_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        lumpedMassTemperatureFvPatchScalarField
    );
}

// ************************************************************************* //
