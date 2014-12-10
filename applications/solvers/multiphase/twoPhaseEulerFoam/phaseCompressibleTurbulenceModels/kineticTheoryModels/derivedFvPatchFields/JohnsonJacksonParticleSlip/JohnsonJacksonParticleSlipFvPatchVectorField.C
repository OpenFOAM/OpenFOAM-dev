/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 OpenFOAM Foundation
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

#include "JohnsonJacksonParticleSlipFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "twoPhaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        JohnsonJacksonParticleSlipFvPatchVectorField
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::JohnsonJacksonParticleSlipFvPatchVectorField::
JohnsonJacksonParticleSlipFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    partialSlipFvPatchVectorField(p, iF),
    specularityCoefficient_(p.size())
{}


Foam::JohnsonJacksonParticleSlipFvPatchVectorField::
JohnsonJacksonParticleSlipFvPatchVectorField
(
    const JohnsonJacksonParticleSlipFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    partialSlipFvPatchVectorField(ptf, p, iF, mapper),
    specularityCoefficient_(ptf.specularityCoefficient_)
{}


Foam::JohnsonJacksonParticleSlipFvPatchVectorField::
JohnsonJacksonParticleSlipFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    partialSlipFvPatchVectorField(p, iF),
    specularityCoefficient_
    (
        "specularityCoefficient",
        dimless,
        dict.lookup("specularityCoefficient")
    )
{
    if
    (
        (specularityCoefficient_.value() < 0)
     || (specularityCoefficient_.value() > 1)
    )
    {
        FatalErrorIn
        (
            "("
                "Foam::JohnsonJacksonParticleSlipFvPatchVectorField::"
                "JohnsonJacksonParticleSlipFvPatchVectorField"
                "const fvPatch& p,"
                "const DimensionedField<scalar, volMesh>& iF,"
                "const dictionary& dict"
            ")"
        )   << "The specularity coefficient has to be between 0 and 1"
            << abort(FatalError);
    }

    fvPatchVectorField::operator=
    (
        vectorField("value", dict, p.size())
    );
}


Foam::JohnsonJacksonParticleSlipFvPatchVectorField::
JohnsonJacksonParticleSlipFvPatchVectorField
(
    const JohnsonJacksonParticleSlipFvPatchVectorField& ptf
)
:
    partialSlipFvPatchVectorField(ptf),
    specularityCoefficient_(ptf.specularityCoefficient_)
{}


Foam::JohnsonJacksonParticleSlipFvPatchVectorField::
JohnsonJacksonParticleSlipFvPatchVectorField
(
    const JohnsonJacksonParticleSlipFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    partialSlipFvPatchVectorField(ptf, iF),
    specularityCoefficient_(ptf.specularityCoefficient_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::JohnsonJacksonParticleSlipFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    partialSlipFvPatchVectorField::autoMap(m);
}


void Foam::JohnsonJacksonParticleSlipFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    partialSlipFvPatchVectorField::rmap(ptf, addr);
}


void Foam::JohnsonJacksonParticleSlipFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // lookup the fluid model and the phase
    const twoPhaseSystem& fluid = db().lookupObject<twoPhaseSystem>
    (
        "phaseProperties"
    );

    const phaseModel& phased
    (
        fluid.phase1().name() == dimensionedInternalField().group()
      ? fluid.phase1()
      : fluid.phase2()
    );

    // lookup all the fields on this patch
    const fvPatchScalarField& alpha
    (
        patch().lookupPatchField<volScalarField, scalar>
        (
            phased.volScalarField::name()
        )
    );

    const fvPatchScalarField& gs0
    (
        patch().lookupPatchField<volScalarField, scalar>
        (
            IOobject::groupName("gs0", phased.name())
        )
    );

    const scalarField nu
    (
        patch().lookupPatchField<volScalarField, scalar>
        (
            IOobject::groupName("nut", phased.name())
        )
    );

    word ThetaName(IOobject::groupName("Theta", phased.name()));

    const fvPatchScalarField& Theta
    (
        db().foundObject<volScalarField>(ThetaName)
      ? patch().lookupPatchField<volScalarField, scalar>(ThetaName)
      : alpha
    );

    // lookup the packed volume fraction
    dimensionedScalar alphaMax
    (
        "alphaMax",
        dimless,
        db()
       .lookupObject<IOdictionary>
        (
            IOobject::groupName("turbulenceProperties", phased.name())
        )
       .subDict("RAS")
       .subDict("kineticTheoryCoeffs")
       .lookup("alphaMax")
    );

    // calculate the slip value fraction
    scalarField c
    (
        constant::mathematical::pi
       *alpha
       *gs0
       *specularityCoefficient_.value()
       *sqrt(3.0*Theta)
       /max(6.0*nu*alphaMax.value(), SMALL)
    );

    this->valueFraction() = c/(c + patch().deltaCoeffs());

    partialSlipFvPatchVectorField::updateCoeffs();
}


void Foam::JohnsonJacksonParticleSlipFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("specularityCoefficient")
        << specularityCoefficient_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// ************************************************************************* //
