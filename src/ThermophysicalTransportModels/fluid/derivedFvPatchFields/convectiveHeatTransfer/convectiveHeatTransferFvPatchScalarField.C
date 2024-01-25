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

#include "convectiveHeatTransferFvPatchScalarField.H"
#include "fluidThermophysicalTransportModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

convectiveHeatTransferFvPatchScalarField::
convectiveHeatTransferFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    L_(dict.lookup<scalar>("L"))
{}


convectiveHeatTransferFvPatchScalarField::
convectiveHeatTransferFvPatchScalarField
(
    const convectiveHeatTransferFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    L_(ptf.L_)
{}


convectiveHeatTransferFvPatchScalarField::
convectiveHeatTransferFvPatchScalarField
(
    const convectiveHeatTransferFvPatchScalarField& htcpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(htcpsf, iF),
    L_(htcpsf.L_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void convectiveHeatTransferFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchi = patch().index();

    const fluidThermophysicalTransportModel& ttm =
        db().lookupType<fluidThermophysicalTransportModel>
        (
            internalField().group()
        );

    const compressibleMomentumTransportModel& turbModel =
        ttm.momentumTransport();

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const scalarField& rhow = turbModel.rho().boundaryField()[patchi];
    const vectorField& Uc = turbModel.U();
    const vectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalarField& Cpw(ttm.thermo().Cp().boundaryField()[patchi]);

    const scalarField kappaEffw(ttm.kappaEff(patchi));
    const scalarField Pr(rhow*nuw*Cpw/kappaEffw);

    scalarField& htc = *this;
    forAll(htc, facei)
    {
        const label celli = patch().faceCells()[facei];

        const scalar Re = mag(Uc[celli] - Uw[facei])*L_/nuw[facei];

        if (Re < 5.0E+05)
        {
            htc[facei] = 0.664*sqrt(Re)*cbrt(Pr[facei])*kappaEffw[facei]/L_;
        }
        else
        {
            htc[facei] = 0.037*pow(Re, 0.8)*cbrt(Pr[facei])*kappaEffw[facei]/L_;
        }
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void convectiveHeatTransferFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeEntry(os, "L", L_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    convectiveHeatTransferFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
