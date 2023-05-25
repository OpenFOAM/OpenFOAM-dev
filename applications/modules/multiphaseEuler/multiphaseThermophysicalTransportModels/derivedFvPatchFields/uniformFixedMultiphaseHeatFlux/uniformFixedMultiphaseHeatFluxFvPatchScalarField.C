/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2023 OpenFOAM Foundation
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

#include "uniformFixedMultiphaseHeatFluxFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "phaseSystem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uniformFixedMultiphaseHeatFluxFvPatchScalarField::
uniformFixedMultiphaseHeatFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    q_(nullptr),
    relax_(1)
{}


Foam::uniformFixedMultiphaseHeatFluxFvPatchScalarField::
uniformFixedMultiphaseHeatFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict, false),
    q_(Function1<scalar>::New("q", dict)),
    relax_(dict.lookupOrDefault<scalar>("relax", 1))
{
    valueFraction() = 1;
    refValue() = patchInternalField();
    refGrad() = Zero;

    operator==(patchInternalField());
}


Foam::uniformFixedMultiphaseHeatFluxFvPatchScalarField::
uniformFixedMultiphaseHeatFluxFvPatchScalarField
(
    const uniformFixedMultiphaseHeatFluxFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(psf, p, iF, mapper),
    q_(psf.q_, false),
    relax_(psf.relax_)
{}


Foam::uniformFixedMultiphaseHeatFluxFvPatchScalarField::
uniformFixedMultiphaseHeatFluxFvPatchScalarField
(
    const uniformFixedMultiphaseHeatFluxFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(psf, iF),
    q_(psf.q_, false),
    relax_(psf.relax_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uniformFixedMultiphaseHeatFluxFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchi = patch().index();

    const scalar q = q_->value(db().time().userTimeValue());

    const phaseSystem& fluid =
        db().lookupObject<phaseSystem>(phaseSystem::propertiesName);

    const phaseModel& thisPhase = fluid.phases()[internalField().group()];

    // Sums of alpha*kappaEff
    scalarField sumAlphaKappaEff(patch().size(), rootVSmall);
    scalarField sumNotThisAlphaKappaEff(patch().size(), rootVSmall);
    scalarField sumNotThisAlphaKappaEffT(patch().size(), rootVSmall);

    // Contributions from phases other than this one
    forAll(fluid.phases(), phasei)
    {
        const phaseModel& phase = fluid.phases()[phasei];

        if (&phase != &thisPhase)
        {
            const scalarField& alpha = phase.boundaryField()[patchi];
            const scalarField kappaEff(phase.kappaEff(patchi));
            const scalarField alphaKappaEff(alpha*kappaEff);
            const fvPatchScalarField& T =
                phase.thermo().T().boundaryField()[patchi];

            sumAlphaKappaEff += alphaKappaEff;
            sumNotThisAlphaKappaEff += alphaKappaEff;
            sumNotThisAlphaKappaEffT += alphaKappaEff*T.patchInternalField();
        }
    }

    // Contribution from this phase, and stabilisation
    const scalarField& alpha = thisPhase.boundaryField()[patchi];
    const scalarField kappaEff(thisPhase.kappaEff(patchi));
    const scalarField alphaKappaEff(alpha*kappaEff);
    const fvPatchScalarField& T =
        thisPhase.thermo().T().boundaryField()[patchi];

    sumAlphaKappaEff += alphaKappaEff;
    sumNotThisAlphaKappaEff =
        max
        (
            sumNotThisAlphaKappaEff,
            rootSmall*kappaEff
        );
    sumNotThisAlphaKappaEffT =
        max
        (
            sumNotThisAlphaKappaEffT,
            rootSmall*kappaEff*T.patchInternalField()
        );

    // Mixed parameters
    valueFraction() = sumNotThisAlphaKappaEff/sumAlphaKappaEff;
    refValue() = sumNotThisAlphaKappaEffT/sumNotThisAlphaKappaEff;
    refGrad() = q/max(alpha, rootSmall)/kappaEff;

    // Modify mixed parameters for under-relaxation
    if (relax_ != 1)
    {
        const scalarField f(valueFraction());
        valueFraction() = 1 - relax_*(1 - f);
        refValue() = (f*relax_*refValue() + (1 - relax_)*T)/valueFraction();
        //refGrad() = refGrad(); // No change
    }

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::uniformFixedMultiphaseHeatFluxFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    writeEntry(os, q_());
    writeEntry(os, "relax", relax_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        uniformFixedMultiphaseHeatFluxFvPatchScalarField
    );
}


// ************************************************************************* //
