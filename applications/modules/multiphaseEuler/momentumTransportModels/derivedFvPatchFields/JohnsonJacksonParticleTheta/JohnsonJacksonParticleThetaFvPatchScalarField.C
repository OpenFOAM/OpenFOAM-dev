/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2024 OpenFOAM Foundation
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

#include "JohnsonJacksonParticleThetaFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "kineticTheoryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        JohnsonJacksonParticleThetaFvPatchScalarField
    );
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::JohnsonJacksonParticleThetaFvPatchScalarField::
JohnsonJacksonParticleThetaFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict, false),
    restitutionCoefficient_
    (
        dict.lookup<scalar>("restitutionCoefficient", unitFraction)
    ),
    specularityCoefficient_
    (
        dict.lookup<scalar>("specularityCoefficient", unitFraction)
    )
{
    if (restitutionCoefficient_ < 0 || restitutionCoefficient_ > 1)
    {
        FatalErrorInFunction
            << "The restitution coefficient has to be between 0 and 1"
            << abort(FatalError);
    }

    if (specularityCoefficient_ < 0 || specularityCoefficient_ > 1)
    {
        FatalErrorInFunction
            << "The specularity coefficient has to be between 0 and 1"
            << abort(FatalError);
    }

    fvPatchScalarField::operator=
    (
        scalarField("value", iF.dimensions(), dict, p.size())
    );
}


Foam::JohnsonJacksonParticleThetaFvPatchScalarField::
JohnsonJacksonParticleThetaFvPatchScalarField
(
    const JohnsonJacksonParticleThetaFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    restitutionCoefficient_(ptf.restitutionCoefficient_),
    specularityCoefficient_(ptf.specularityCoefficient_)
{
}


Foam::JohnsonJacksonParticleThetaFvPatchScalarField::
JohnsonJacksonParticleThetaFvPatchScalarField
(
    const JohnsonJacksonParticleThetaFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    restitutionCoefficient_(ptf.restitutionCoefficient_),
    specularityCoefficient_(ptf.specularityCoefficient_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::JohnsonJacksonParticleThetaFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // lookup the fluid model and the phase
    const phaseSystem& fluid =
        db().lookupObject<phaseSystem>(phaseSystem::propertiesName);

    const phaseModel& phase
    (
        fluid.phases()[internalField().group()]
    );

    // lookup all the fields on this patch
    const fvPatchScalarField& alpha
    (
        patch().lookupPatchField<volScalarField, scalar>
        (
            phase.volScalarField::name()
        )
    );

    const fvPatchVectorField& U
    (
        patch().lookupPatchField<volVectorField, vector>
        (
            IOobject::groupName("U", phase.name())
        )
    );

    const fvPatchScalarField& gs0
    (
        patch().lookupPatchField<volScalarField, scalar>
        (
            IOobject::groupName
            (
                Foam::typedName<RASModels::kineticTheoryModel>("gs0"),
                phase.name()
            )
        )
    );

    const fvPatchScalarField& kappa
    (
        patch().lookupPatchField<volScalarField, scalar>
        (
            IOobject::groupName
            (
                Foam::typedName<RASModels::kineticTheoryModel>("kappa"),
                phase.name()
            )
        )
    );

    const scalarField Theta(patchInternalField());

    // calculate the reference value and the value fraction
    if (restitutionCoefficient_ != 1.0)
    {
        this->refValue() =
            (2.0/3.0)
           *specularityCoefficient_
           *magSqr(U)
           /(scalar(1) - sqr(restitutionCoefficient_));

        this->refGrad() = 0.0;

        scalarField c
        (
             constant::mathematical::pi
            *alpha
            *gs0
            *(scalar(1) - sqr(restitutionCoefficient_))
            *sqrt(3*Theta)
            /max(4*kappa*phase.alphaMax(), small)
        );

        this->valueFraction() = c/(c + patch().deltaCoeffs());
    }

    // for a restitution coefficient of 1, the boundary degenerates to a fixed
    // gradient condition
    else
    {
        this->refValue() = 0.0;

        this->refGrad() =
            pos0(alpha - small)
           *constant::mathematical::pi
           *specularityCoefficient_
           *alpha
           *gs0
           *sqrt(3*Theta)
           *magSqr(U)
           /max(6*kappa*phase.alphaMax(), small);

        this->valueFraction() = 0;
    }

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::JohnsonJacksonParticleThetaFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    writeEntry(os, "restitutionCoefficient", restitutionCoefficient_);
    writeEntry(os, "specularityCoefficient", specularityCoefficient_);
    writeEntry(os, "value", *this);
}


// ************************************************************************* //
