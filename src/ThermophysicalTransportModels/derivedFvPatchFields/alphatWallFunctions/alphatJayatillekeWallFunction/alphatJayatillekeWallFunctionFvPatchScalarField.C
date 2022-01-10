/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
#include "thermophysicalTransportModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

scalar alphatJayatillekeWallFunctionFvPatchScalarField::maxExp_ = 50.0;
scalar alphatJayatillekeWallFunctionFvPatchScalarField::tolerance_ = 0.01;
label alphatJayatillekeWallFunctionFvPatchScalarField::maxIters_ = 10;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

scalar alphatJayatillekeWallFunctionFvPatchScalarField::Psmooth
(
    const scalar Prat
) const
{
    return 9.24*(pow(Prat, 0.75) - 1.0)*(1.0 + 0.28*exp(-0.007*Prat));
}


scalar alphatJayatillekeWallFunctionFvPatchScalarField::yPlusTherm
(
    const nutWallFunctionFvPatchScalarField& nutw,
    const scalar P,
    const scalar Prat
) const
{
    scalar ypt = 11.0;

    for (int i=0; i<maxIters_; i++)
    {
        scalar f = ypt - (log(nutw.E()*ypt)/nutw.kappa() + P)/Prat;
        scalar df = 1.0 - 1.0/(ypt*nutw.kappa()*Prat);
        scalar yptNew = ypt - f/df;

        if (yptNew < vSmall)
        {
            return 0;
        }
        else if (mag(yptNew - ypt) < tolerance_)
        {
            return yptNew;
        }
        else
        {
            ypt = yptNew;
        }
     }

    return ypt;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

alphatJayatillekeWallFunctionFvPatchScalarField::
alphatJayatillekeWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    Prt_(0.85)
{}


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
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    Prt_(ptf.Prt_)
{}


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

void alphatJayatillekeWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchi = patch().index();

    const thermophysicalTransportModel& ttm =
        db().lookupObject<thermophysicalTransportModel>
        (
            IOobject::groupName
            (
                thermophysicalTransportModel::typeName,
                internalField().group()
            )
        );

    const compressibleMomentumTransportModel& turbModel =
        ttm.momentumTransport();

    const nutWallFunctionFvPatchScalarField& nutw =
        nutWallFunctionFvPatchScalarField::nutw(turbModel, patchi);

    const scalar Cmu25 = pow025(nutw.Cmu());

    const scalarField& y = turbModel.y()[patchi];

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const tmp<scalarField> talphaw
    (
        ttm.thermo().kappa().boundaryField()[patchi]
       /ttm.thermo().Cp().boundaryField()[patchi]
    );
    const scalarField& alphaw = talphaw();

    scalarField& alphatw = *this;

    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalarField magUp(mag(Uw.patchInternalField() - Uw));
    const scalarField magGradUw(mag(Uw.snGrad()));

    const scalarField& rhow = turbModel.rho().boundaryField()[patchi];
    const fvPatchScalarField& hew = ttm.thermo().he().boundaryField()[patchi];

    // Heat flux [W/m^2] - lagging alphatw
    const scalarField qDot(ttm.thermo().alphaEff(alphatw, patchi)*hew.snGrad());

    // Populate boundary values
    forAll(alphatw, facei)
    {
        label celli = patch().faceCells()[facei];

        scalar uTau = Cmu25*sqrt(k[celli]);

        scalar yPlus = uTau*y[facei]/nuw[facei];

        // Molecular Prandtl number
        scalar Pr = rhow[facei]*nuw[facei]/alphaw[facei];

        // Molecular-to-turbulent Prandtl number ratio
        scalar Prat = Pr/Prt_;

        // Thermal sublayer thickness
        scalar P = Psmooth(Prat);
        scalar yPlusTherm = this->yPlusTherm(nutw, P, Prat);

        // Evaluate new effective thermal diffusivity
        scalar alphaEff = 0.0;
        if (yPlus < yPlusTherm)
        {
            const scalar A = qDot[facei]*rhow[facei]*uTau*y[facei];

            const scalar B = qDot[facei]*Pr*yPlus;

            const scalar C = Pr*0.5*rhow[facei]*uTau*sqr(magUp[facei]);

            alphaEff = A/(B + C + vSmall);
        }
        else
        {
            const scalar A = qDot[facei]*rhow[facei]*uTau*y[facei];

            const scalar B =
                qDot[facei]*Prt_*(1.0/nutw.kappa()*log(nutw.E()*yPlus) + P);

            const scalar magUc =
                uTau/nutw.kappa()*log(nutw.E()*yPlusTherm) - mag(Uw[facei]);

            const scalar C =
                0.5*rhow[facei]*uTau
               *(Prt_*sqr(magUp[facei]) + (Pr - Prt_)*sqr(magUc));

            alphaEff = A/(B + C + vSmall);
        }

        // Update turbulent thermal diffusivity
        alphatw[facei] = max(0.0, alphaEff - alphaw[facei]);

        if (debug)
        {
            Info<< "    uTau           = " << uTau << nl
                << "    Pr             = " << Pr << nl
                << "    Prt            = " << Prt_ << nl
                << "    qDot           = " << qDot[facei] << nl
                << "    yPlus          = " << yPlus << nl
                << "    yPlusTherm     = " << yPlusTherm << nl
                << "    alphaEff       = " << alphaEff << nl
                << "    alphaw         = " << alphaw[facei] << nl
                << "    alphatw        = " << alphatw[facei] << nl
                << endl;
        }
    }

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
