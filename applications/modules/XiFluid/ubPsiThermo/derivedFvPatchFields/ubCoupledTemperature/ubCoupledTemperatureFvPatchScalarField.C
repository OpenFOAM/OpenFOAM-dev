/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2025 OpenFOAM Foundation
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

#include "ubCoupledTemperatureFvPatchScalarField.H"
#include "XiFluid.H"
#include "fieldMapper.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::ubCoupledTemperatureFvPatchScalarField::getThis
(
    tmp<scalarField>& kappa,
    tmp<scalarField>& sumKappaTByDelta,
    tmp<scalarField>& sumKappaByDelta,
    scalarField& sumq,
    tmp<scalarField>& qByKappa
) const
{
    const solvers::XiFluid& XiFluid
    (
        patch().boundaryMesh().mesh()
       .lookupObject<solvers::XiFluid>(solver::typeName)
    );

    const ubPsiThermo& thermo = XiFluid.thermo;
    const ubPsiMulticomponentThermo& uThermo = thermo.uThermo();
    const ubPsiMulticomponentThermo& bThermo = thermo.bThermo();

    const scalarField b =
          thermo.rho(patch().index())
         *thermo.b().boundaryField()[patch().index()]
         /uThermo.rho(patch().index());

    const scalarField c =
          thermo.rho(patch().index())
         *thermo.c().boundaryField()[patch().index()]
         /bThermo.rho(patch().index());

    const fvPatchScalarField& Twu(uThermo.T().boundaryField()[patch().index()]);

    const fvPatchScalarField& Twb(bThermo.T().boundaryField()[patch().index()]);

    const scalarField nut(XiFluid.momentumTransport.nut(patch().index()));

    // Simple model pending generalised thermophysicalTransportModel interface
    tmp<scalarField> kappaEffu
    (
        uThermo.rho(patch().index())*nut
      + uThermo.mu().boundaryField()[patch().index()]
    );
    tmp<scalarField> kappaEffb
    (
        bThermo.rho(patch().index())*nut
      + bThermo.mu().boundaryField()[patch().index()]
    );

    tmp<scalarField> alphaKappaEffu(b*kappaEffu());
    tmp<scalarField> alphaKappaEffb(c*kappaEffb());

    scalarField sumKappa(size(), scalar(0));
    scalarField sumKappaT(size(), scalar(0));

    if (&Twu == this)
    {
        kappa = alphaKappaEffu;
        qByKappa = sumq/kappaEffu;
        sumq -= b*sumq;
    }
    else
    {
        const scalarField Tu
        (
            uThermo.T().boundaryField()[patch().index()]
           .patchInternalField()
        );

        sumKappa += alphaKappaEffu();
        sumKappaT += alphaKappaEffu*Tu;
    }

    if (&Twb == this)
    {
        kappa = alphaKappaEffb;
        qByKappa = sumq/kappaEffb;
        sumq -= c*sumq;
    }
    else
    {
        const scalarField Tb
        (
            bThermo.T().boundaryField()[patch().index()]
           .patchInternalField()
        );

        sumKappa += alphaKappaEffb();
        sumKappaT += alphaKappaEffb*Tb;
    }

    sumKappaByDelta = sumKappa*patch().deltaCoeffs();
    sumKappaTByDelta = sumKappaT*patch().deltaCoeffs();
}


void Foam::ubCoupledTemperatureFvPatchScalarField::getNbr
(
    tmp<scalarField>& sumKappaTByDeltaNbr,
    tmp<scalarField>& sumKappaByDeltaNbr,
    tmp<scalarField>& qNbr
) const
{
    const solvers::XiFluid& XiFluid
    (
        patch().boundaryMesh().mesh()
       .lookupObject<solvers::XiFluid>(solver::typeName)
    );

    const ubPsiThermo& thermo = XiFluid.thermo;
    const ubPsiMulticomponentThermo& uThermo = thermo.uThermo();
    const ubPsiMulticomponentThermo& bThermo = thermo.bThermo();

    const scalarField b =
          thermo.rho(patch().index())
         *thermo.b().boundaryField()[patch().index()]
         /uThermo.rho(patch().index());

    const scalarField c =
          thermo.rho(patch().index())
         *thermo.c().boundaryField()[patch().index()]
         /bThermo.rho(patch().index());

    const scalarField Tu
    (
        uThermo.T().boundaryField()[patch().index()].patchInternalField()
    );

    const scalarField Tb
    (
        bThermo.T().boundaryField()[patch().index()].patchInternalField()
    );

    const scalarField nut(XiFluid.momentumTransport.nut(patch().index()));

    // Simple model pending generalised thermophysicalTransportModel interface
    tmp<scalarField> kappaEffu
    (
        uThermo.rho(patch().index())*nut
      + uThermo.mu().boundaryField()[patch().index()]
    );
    tmp<scalarField> kappaEffb
    (
        bThermo.rho(patch().index())*nut
      + bThermo.mu().boundaryField()[patch().index()]
    );

    const scalarField alphaKappaEffu(b*kappaEffu);
    const scalarField alphaKappaEffb(c*kappaEffb);

    sumKappaByDeltaNbr =
        (alphaKappaEffu + alphaKappaEffb)*patch().deltaCoeffs();
    sumKappaTByDeltaNbr =
        (alphaKappaEffu*Tu + alphaKappaEffb*Tb)*patch().deltaCoeffs();
}


void Foam::ubCoupledTemperatureFvPatchScalarField::getNbr
(
    tmp<scalarField>& TrefNbr,
    tmp<scalarField>& qNbr
) const
{
    const solvers::XiFluid& XiFluid
    (
        patch().boundaryMesh().mesh()
       .lookupObject<solvers::XiFluid>(solver::typeName)
    );

    const ubPsiThermo& thermo = XiFluid.thermo;
    const ubPsiMulticomponentThermo& uThermo = thermo.uThermo();
    const ubPsiMulticomponentThermo& bThermo = thermo.bThermo();

    const scalarField b =
          thermo.rho(patch().index())
         *thermo.b().boundaryField()[patch().index()]
         /uThermo.rho(patch().index());

    const scalarField c =
          thermo.rho(patch().index())
         *thermo.c().boundaryField()[patch().index()]
         /bThermo.rho(patch().index());

    TrefNbr =
       b*uThermo.T().boundaryField()[patch().index()]
     + c*bThermo.T().boundaryField()[patch().index()];
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ubCoupledTemperatureFvPatchScalarField::
ubCoupledTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    coupledTemperatureFvPatchScalarField(p, iF, dict)
{}


Foam::ubCoupledTemperatureFvPatchScalarField::
ubCoupledTemperatureFvPatchScalarField
(
    const ubCoupledTemperatureFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    coupledTemperatureFvPatchScalarField(psf, p, iF, mapper)
{}


Foam::ubCoupledTemperatureFvPatchScalarField::
ubCoupledTemperatureFvPatchScalarField
(
    const ubCoupledTemperatureFvPatchScalarField& psf,
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
        ubCoupledTemperatureFvPatchScalarField
    );
}


// ************************************************************************* //
