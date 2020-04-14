/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
#include "thermophysicalTransportModel.H"
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
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    L_(1.0)
{}


convectiveHeatTransferFvPatchScalarField::
convectiveHeatTransferFvPatchScalarField
(
    const convectiveHeatTransferFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    L_(ptf.L_)
{}


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
    const convectiveHeatTransferFvPatchScalarField& htcpsf
)
:
    fixedValueFvPatchScalarField(htcpsf),
    L_(htcpsf.L_)
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

    const scalarField alphaEffw(ttm.alphaEff(patchi));

    const tmp<scalarField> tmuw = turbModel.mu(patchi);
    const scalarField& muw = tmuw();

    const scalarField& rhow = turbModel.rho().boundaryField()[patchi];
    const vectorField& Uc = turbModel.U();
    const vectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalarField& Tw = ttm.thermo().T().boundaryField()[patchi];
    const scalarField Cpw(ttm.thermo().Cp(Tw, patchi));

    const scalarField kappaw(Cpw*alphaEffw);
    const scalarField Pr(muw*Cpw/kappaw);

    scalarField& htc = *this;
    forAll(htc, facei)
    {
        label celli = patch().faceCells()[facei];

        scalar Re = rhow[facei]*mag(Uc[celli] - Uw[facei])*L_/muw[facei];

        if (Re < 5.0E+05)
        {
            htc[facei] = 0.664*sqrt(Re)*cbrt(Pr[facei])*kappaw[facei]/L_;
        }
        else
        {
            htc[facei] = 0.037*pow(Re, 0.8)*cbrt(Pr[facei])*kappaw[facei]/L_;
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
