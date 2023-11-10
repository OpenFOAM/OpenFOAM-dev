/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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

#include "coupledMultiphaseTemperatureFvPatchScalarField.H"
#include "fieldMapper.H"
#include "phaseSystem.H"
#include "compressibleMomentumTransportModel.H"
#include "phaseCompressibleMomentumTransportModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::coupledMultiphaseTemperatureFvPatchScalarField::getThis
(
    tmp<scalarField>& kappa,
    tmp<scalarField>& sumKappaTByDelta,
    tmp<scalarField>& sumKappaByDelta,
    scalarField& sumq,
    tmp<scalarField>& qByKappa
) const
{
    // Lookup the fluid model
    const phaseSystem& fluid =
        patch().boundaryMesh().mesh()
       .lookupObject<phaseSystem>(phaseSystem::propertiesName);

    scalarField sumKappa(size(), scalar(0));
    scalarField sumKappaT(size(), scalar(0));

    forAll(fluid.phases(), phasei)
    {
        const phaseModel& phase = fluid.phases()[phasei];
        const rhoThermo& thermo = phase.thermo();

        const fvPatchScalarField& Tw =
            thermo.T().boundaryField()[patch().index()];

        const fvPatchScalarField& alpha =
            phase.boundaryField()[patch().index()];

        tmp<scalarField> kappaEff(phase.kappaEff(patch().index()));
        tmp<scalarField> alphaKappaEff(alpha*kappaEff());

        if (&Tw == this)
        {
            kappa = alphaKappaEff;
            qByKappa = sumq/kappaEff;
            sumq -= alpha*sumq;
        }
        else
        {
            const scalarField T
            (
                thermo.T().boundaryField()[patch().index()]
               .patchInternalField()
            );

            sumKappa += alphaKappaEff();
            sumKappaT += alphaKappaEff*T;
        }
    }

    sumKappaByDelta = sumKappa*patch().deltaCoeffs();
    sumKappaTByDelta = sumKappaT*patch().deltaCoeffs();
}


void Foam::coupledMultiphaseTemperatureFvPatchScalarField::getNbr
(
    tmp<scalarField>& sumKappaTByDeltaNbr,
    tmp<scalarField>& sumKappaByDeltaNbr,
    tmp<scalarField>& qNbr
) const
{
    // Lookup the fluid model
    const phaseSystem& fluid =
        patch().boundaryMesh().mesh()
       .lookupObject<phaseSystem>(phaseSystem::propertiesName);

    scalarField sumKappa(size(), scalar(0));
    scalarField sumKappaT(size(), scalar(0));

    forAll(fluid.phases(), phasei)
    {
        const phaseModel& phase = fluid.phases()[phasei];
        const rhoThermo& thermo = phase.thermo();

        const fvPatchScalarField& alpha =
            phase.boundaryField()[patch().index()];

        const scalarField T
        (
            thermo.T().boundaryField()[patch().index()].patchInternalField()
        );

        const scalarField alphaKappaEff(alpha*phase.kappaEff(patch().index()));

        sumKappa += alphaKappaEff;
        sumKappaT += alphaKappaEff*T;
    }

    sumKappaByDeltaNbr = sumKappa*patch().deltaCoeffs();
    sumKappaTByDeltaNbr = sumKappaT*patch().deltaCoeffs();
}


void Foam::coupledMultiphaseTemperatureFvPatchScalarField::getNbr
(
    tmp<scalarField>& TrefNbr,
    tmp<scalarField>& qNbr
) const
{
    // Lookup the fluid model
    const phaseSystem& fluid =
        patch().boundaryMesh().mesh()
       .lookupObject<phaseSystem>(phaseSystem::propertiesName);

    TrefNbr = new scalarField(size(), scalar(0));
    scalarField& Tw = TrefNbr.ref();

    forAll(fluid.phases(), phasei)
    {
        const phaseModel& phase = fluid.phases()[phasei];
        const rhoThermo& thermo = phase.thermo();

        const fvPatchScalarField& alpha =
            phase.boundaryField()[patch().index()];

        Tw += alpha*thermo.T().boundaryField()[patch().index()];

        /*
        // Pending the addition of phase anisotropic thermal transport
        tmp<scalarField> qCorr = ttm.qCorr(patch().index());

        if (qCorr.valid())
        {
            if (qCorr.valid())
            {
                qNbr.ref() += alpha*qCorr;
            }
            else
            {
                qNbr = alpha*qCorr;
            }
        }
        */
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coupledMultiphaseTemperatureFvPatchScalarField::
coupledMultiphaseTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    coupledTemperatureFvPatchScalarField(p, iF, dict)
{}


Foam::coupledMultiphaseTemperatureFvPatchScalarField::
coupledMultiphaseTemperatureFvPatchScalarField
(
    const coupledMultiphaseTemperatureFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    coupledTemperatureFvPatchScalarField(psf, p, iF, mapper)
{}


Foam::coupledMultiphaseTemperatureFvPatchScalarField::
coupledMultiphaseTemperatureFvPatchScalarField
(
    const coupledMultiphaseTemperatureFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    coupledTemperatureFvPatchScalarField(psf, iF)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        coupledMultiphaseTemperatureFvPatchScalarField
    );
}


// ************************************************************************* //
