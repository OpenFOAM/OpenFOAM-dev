/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2025 OpenFOAM Foundation
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

#include "multiphaseExternalTemperatureFvPatchScalarField.H"
#include "phaseSystem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::multiphaseExternalTemperatureFvPatchScalarField::getKappa
(
    scalarField& kappa,
    tmp<scalarField>& sumKappaTcByDelta,
    tmp<scalarField>& sumKappaByDelta,
    tmp<scalarField>& T,
    tmp<scalarField>& sumq
) const
{
    const phaseSystem& fluid =
        db().lookupObject<phaseSystem>(phaseSystem::propertiesName);

    const label patchi = patch().index();

    const phaseModel& thisPhase = fluid.phases()[internalField().group()];

    // Get alpha and kappa for this phase
    const scalarField& alpha = thisPhase.boundaryField()[patchi];
    const scalarField alphaKappaEff(alpha*thisPhase.kappaEff(patchi));

    // Get alpha and kappa sums for all other phases
    scalarField sumKappa(size(), scalar(0));
    scalarField sumKappaTc(size(), scalar(0));
    scalarField sumKappaT(size(), scalar(0));
    forAll(fluid.phases(), phasei)
    {
        const phaseModel& phase = fluid.phases()[phasei];

        if (&phase != &thisPhase)
        {
            const scalarField& alpha = phase.boundaryField()[patchi];
            const scalarField alphaKappaEff(alpha*phase.kappaEff(patchi));

            const fvPatchScalarField& T =
                phase.thermo().T().boundaryField()[patchi];

            sumKappa += alphaKappaEff;
            sumKappaTc += alphaKappaEff*T.patchInternalField();
            sumKappaT += alphaKappaEff*T;
        }
    }

    kappa = alphaKappaEff;

    if (fluid.phases().size() > 1)
    {
        sumKappaByDelta = sumKappa*patch().deltaCoeffs();
        sumKappaTcByDelta = sumKappaTc*patch().deltaCoeffs();
    }

    T = (sumKappaT + alphaKappaEff*(*this))/(sumKappa + alphaKappaEff);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiphaseExternalTemperatureFvPatchScalarField::
multiphaseExternalTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    externalTemperatureFvPatchScalarField(p, iF, dict)
{}


Foam::multiphaseExternalTemperatureFvPatchScalarField::
multiphaseExternalTemperatureFvPatchScalarField
(
    const multiphaseExternalTemperatureFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    externalTemperatureFvPatchScalarField(psf, p, iF, mapper)
{}


Foam::multiphaseExternalTemperatureFvPatchScalarField::
multiphaseExternalTemperatureFvPatchScalarField
(
    const multiphaseExternalTemperatureFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    externalTemperatureFvPatchScalarField(psf, iF)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        multiphaseExternalTemperatureFvPatchScalarField
    );
}


// ************************************************************************* //
