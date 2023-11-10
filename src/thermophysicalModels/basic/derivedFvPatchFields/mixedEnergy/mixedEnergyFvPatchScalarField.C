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

#include "mixedEnergyFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fieldMapper.H"
#include "volFields.H"
#include "basicThermo.H"
#include "mixedEnergyCalculatedTemperatureFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixedEnergyFvPatchScalarField::
mixedEnergyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF)
{
    valueFraction() = 0.0;
    refValue() = 0.0;
    refGrad() = 0.0;
}


Foam::mixedEnergyFvPatchScalarField::
mixedEnergyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict)
{}


Foam::mixedEnergyFvPatchScalarField::
mixedEnergyFvPatchScalarField
(
    const mixedEnergyFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::mixedEnergyFvPatchScalarField::
mixedEnergyFvPatchScalarField
(
    const mixedEnergyFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mixedEnergyFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const basicThermo& thermo = basicThermo::lookupThermo(*this);
    const label patchi = patch().index();

    fvPatchScalarField& Tp =
        const_cast<fvPatchScalarField&>(thermo.T().boundaryField()[patchi]);

    if (isA<mixedFvPatchScalarField>(Tp))
    {
        mixedFvPatchScalarField& Tm =
            refCast<mixedFvPatchScalarField>(Tp);

        Tm.evaluate();

        valueFraction() = Tm.valueFraction();
        refValue() = thermo.he(Tm.refValue(), patchi);
        refGrad() =
            thermo.Cpv(Tm, patchi)*Tm.refGrad()
          + patch().deltaCoeffs()*
            (
                thermo.he(Tm, patchi)
              - thermo.he(Tm, patch().faceCells())
            );
    }
    else if (isA<mixedEnergyCalculatedTemperatureFvPatchScalarField>(Tp))
    {
        mixedEnergyCalculatedTemperatureFvPatchScalarField& Tm =
            refCast<mixedEnergyCalculatedTemperatureFvPatchScalarField>(Tp);

        Tm.evaluate();

        valueFraction() = Tm.heValueFraction();
        refValue() = Tm.heRefValue();
        refGrad() = Tm.heRefGrad();
    }
    else
    {
        FatalErrorInFunction
            << "Temperature boundary condition not recognised."
            << "A " << typeName << " condition for energy must be used with a "
            << mixedFvPatchScalarField::typeName << " or "
            << mixedEnergyCalculatedTemperatureFvPatchScalarField::typeName
            << " condition for temperature."
            << exit(FatalError);
    }

    mixedFvPatchScalarField::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeNullConstructablePatchTypeField
    (
        fvPatchScalarField,
        mixedEnergyFvPatchScalarField
    );
}

// ************************************************************************* //
