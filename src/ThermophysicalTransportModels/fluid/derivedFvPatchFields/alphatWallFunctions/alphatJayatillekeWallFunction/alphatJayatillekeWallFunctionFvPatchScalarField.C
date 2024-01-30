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

#include "alphatJayatillekeWallFunctionFvPatchScalarField.H"
#include "fluidThermophysicalTransportModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * Private Static Data Members * * * * * * * * * * * //

const scalar alphatJayatillekeWallFunctionFvPatchScalarField::tolerance_ = 0.01;
const label alphatJayatillekeWallFunctionFvPatchScalarField::maxIters_ = 10;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

alphatJayatillekeWallFunctionFvPatchScalarField::
alphatJayatillekeWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    Prt_(dict.lookupOrDefault<scalar>("Prt", 0.85))
{}


alphatJayatillekeWallFunctionFvPatchScalarField::
alphatJayatillekeWallFunctionFvPatchScalarField
(
    const alphatJayatillekeWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper, false),
    Prt_(ptf.Prt_)
{
    mapper(*this, ptf, [&](){ return this->patchInternalField(); });
}


alphatJayatillekeWallFunctionFvPatchScalarField::
alphatJayatillekeWallFunctionFvPatchScalarField
(
    const alphatJayatillekeWallFunctionFvPatchScalarField& awfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(awfpsf, iF),
    Prt_(awfpsf.Prt_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void alphatJayatillekeWallFunctionFvPatchScalarField::map
(
    const fvPatchScalarField& ptf,
    const fieldMapper& mapper
)
{
    mapper(*this, ptf, [&](){ return this->patchInternalField(); });
}


tmp<scalarField> alphatJayatillekeWallFunctionFvPatchScalarField::P
(
    const scalarField& Prat
)
{
    return 9.24*(pow(Prat, 0.75) - 1)*(1 + 0.28*exp(-0.007*Prat));
}


tmp<scalarField> alphatJayatillekeWallFunctionFvPatchScalarField::yPlusTherm
(
    const nutWallFunctionFvPatchScalarField& nutw,
    const scalarField& P,
    const scalarField& Prat
)
{
    tmp<scalarField> typt(new scalarField(nutw.size()));
    scalarField& ypt = typt.ref();

    const scalar E = nutw.E();
    const scalar kappa = nutw.kappa();

    forAll(ypt, facei)
    {
        ypt[facei] = 11;

        for (int i=0; i<maxIters_; i++)
        {
            const scalar f =
                ypt[facei] - (log(E*ypt[facei])/kappa + P[facei])/Prat[facei];

            const scalar df = 1 - 1/(ypt[facei]*nutw.kappa()*Prat[facei]);

            const scalar dypt = - f/df;

            ypt[facei] += dypt;

            if (ypt[facei] < vSmall)
            {
                ypt[facei] = 0;
                break;
            }

            if (mag(dypt) < tolerance_)
            {
                break;
            }
        }
    }

    return typt;
}


tmp<scalarField> alphatJayatillekeWallFunctionFvPatchScalarField::alphat
(
    const fluidThermophysicalTransportModel& ttm,
    const scalar Prt,
    const label patchi
)
{
    const compressibleMomentumTransportModel& turbModel =
        ttm.momentumTransport();

    const nutWallFunctionFvPatchScalarField& nutw =
        nutWallFunctionFvPatchScalarField::nutw(turbModel, patchi);

    const scalar Cmu25 = pow025(nutw.Cmu());

    const scalarField& y = turbModel.yb()[patchi];

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const scalarField& Cpw(ttm.thermo().Cp().boundaryField()[patchi]);
    const scalarField alphaw(ttm.thermo().kappa().boundaryField()[patchi]/Cpw);

    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalarField magUp(mag(Uw.patchInternalField() - Uw));
    const scalarField magGradUw(mag(Uw.snGrad()));

    const scalarField& rhow = turbModel.rho().boundaryField()[patchi];
    const fvPatchScalarField& Tw = ttm.thermo().T().boundaryField()[patchi];

    // Temperature gradient
    const scalarField gradTw(Tw.snGrad());

    // Molecular Prandtl number
    const scalarField Pr(rhow*nuw/alphaw);

    // Molecular-to-turbulent Prandtl number ratio
    const scalarField Prat(Pr/Prt);

    // Thermal sublayer thickness
    const scalarField P
    (
        alphatJayatillekeWallFunctionFvPatchScalarField::P(Prat)
    );
    const scalarField yPlusTherm
    (
        alphatJayatillekeWallFunctionFvPatchScalarField::yPlusTherm
        (
            nutw,
            P,
            Prat
        )
    );

    // Populate boundary values
    tmp<scalarField> talphatw(new scalarField(nutw.size()));
    scalarField& alphatw = talphatw.ref();
    forAll(alphatw, facei)
    {
        const label celli = nutw.patch().faceCells()[facei];
        const scalar uTau = Cmu25*sqrt(k[celli]);
        const scalar yPlus = uTau*y[facei]/nuw[facei];

        // Evaluate new effective thermal diffusivity
        scalar alphaEff = 0;
        if (yPlus < yPlusTherm[facei])
        {
            const scalar A = Cpw[facei]*gradTw[facei]*rhow[facei]*uTau*y[facei];
            const scalar B = Cpw[facei]*gradTw[facei]*Pr[facei]*yPlus;
            const scalar C = Pr[facei]*0.5*rhow[facei]*uTau*sqr(magUp[facei]);

            alphaEff = (A - C)/(B + sign(B)*rootVSmall);
        }
        else
        {
            const scalar A = Cpw[facei]*gradTw[facei]*rhow[facei]*uTau*y[facei];

            const scalar B =
                Cpw[facei]*gradTw[facei]*Prt
               *(1/nutw.kappa()*log(nutw.E()*yPlus) + P[facei]);

            const scalar magUc =
                uTau/nutw.kappa()
               *log(nutw.E()*yPlusTherm[facei]) - mag(Uw[facei]);

            const scalar C =
                0.5*rhow[facei]*uTau
               *(Prt*sqr(magUp[facei]) + (Pr[facei] - Prt)*sqr(magUc));

            alphaEff = (A - C)/(B + sign(B)*rootVSmall);
        }

        // Bounds on turbulent thermal diffusivity
        static const scalar alphatwMin = 0;
        const scalar alphatwMax = great*alphaw[facei]*nutw[facei]/nuw[facei];

        // Update turbulent thermal diffusivity
        alphatw[facei] =
            min(max(alphaEff - alphaw[facei], alphatwMin), alphatwMax);
    }

    return talphatw;
}


void alphatJayatillekeWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fluidThermophysicalTransportModel& ttm =
        db().lookupType<fluidThermophysicalTransportModel>
        (
            internalField().group()
        );

    this->operator==(alphat(ttm, Prt_, patch().index()));

    fixedValueFvPatchField<scalar>::updateCoeffs();
}


void alphatJayatillekeWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeEntry(os, "Prt", Prt_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    alphatJayatillekeWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
