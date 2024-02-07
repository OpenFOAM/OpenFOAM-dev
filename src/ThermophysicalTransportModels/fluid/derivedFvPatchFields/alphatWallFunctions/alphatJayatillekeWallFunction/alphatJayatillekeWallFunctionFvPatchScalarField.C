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

    const scalar E = nutw.E();
    const scalar kappa = nutw.kappa();

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const scalarField alphaw
    (
        ttm.thermo().kappa().boundaryField()[patchi]
       /ttm.thermo().Cp().boundaryField()[patchi]
    );

    const scalarField& rhow = turbModel.rho().boundaryField()[patchi];

    // Molecular Prandtl number
    const scalarField Pr(rhow*nuw/alphaw);

    // Molecular-to-turbulent Prandtl number ratio
    const scalarField Prat(Pr/Prt);

    // Momentum sublayer thickness
    const scalarField yPlus(nutw.yPlus());

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
    tmp<scalarField> talphatw(new scalarField(nutw.size(), scalar(0)));
    scalarField& alphatw = talphatw.ref();
    forAll(alphatw, facei)
    {
        if (yPlus[facei] > yPlusTherm[facei])
        {
            const scalar Tplus = Prt*(log(E*yPlus[facei])/kappa + P[facei]);

            const scalar alphaByAlphaEff = Tplus/Pr[facei]/yPlus[facei];

            alphatw[facei] =
                alphaw[facei]*max(1/alphaByAlphaEff - 1, scalar(0));
        }
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
