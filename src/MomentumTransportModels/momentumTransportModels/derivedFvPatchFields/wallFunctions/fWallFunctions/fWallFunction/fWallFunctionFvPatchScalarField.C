/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2020 OpenFOAM Foundation
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

#include "fWallFunctionFvPatchScalarField.H"
#include "nutWallFunctionFvPatchScalarField.H"
#include "momentumTransportModel.H"
#include "v2f.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fWallFunctionFvPatchScalarField::fWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF)
{}


fWallFunctionFvPatchScalarField::fWallFunctionFvPatchScalarField
(
    const fWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper)
{}


fWallFunctionFvPatchScalarField::fWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict)
{}


fWallFunctionFvPatchScalarField::fWallFunctionFvPatchScalarField
(
    const fWallFunctionFvPatchScalarField& v2wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(v2wfpsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void fWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchi = patch().index();

    const momentumTransportModel& turbModel =
        db().lookupObject<momentumTransportModel>
        (
            IOobject::groupName
            (
                momentumTransportModel::typeName,
                internalField().group()
            )
        );
    const v2fBase& v2fModel = refCast<const v2fBase>(turbModel);

    const nutWallFunctionFvPatchScalarField& nutw =
        nutWallFunctionFvPatchScalarField::nutw(turbModel, patchi);

    const scalarField& y = turbModel.y()[patchi];

    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();

    const tmp<volScalarField> tepsilon = turbModel.epsilon();
    const volScalarField& epsilon = tepsilon();

    const tmp<volScalarField> tv2 = v2fModel.v2();
    const volScalarField& v2 = tv2();

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const scalar Cmu25 = pow025(nutw.Cmu());

    scalarField& f = *this;

    // Set f wall values
    forAll(f, facei)
    {
        label celli = patch().faceCells()[facei];

        scalar uTau = Cmu25*sqrt(k[celli]);

        scalar yPlus = uTau*y[facei]/nuw[facei];

        if (yPlus > nutw.yPlusLam())
        {
            scalar N = 6.0;
            scalar v2c = v2[celli];
            scalar epsc = epsilon[celli];
            scalar kc = k[celli];

            f[facei] = N*v2c*epsc/(sqr(kc) + rootVSmall);
            f[facei] /= sqr(uTau) + rootVSmall;
        }
        else
        {
            f[facei] = 0.0;
        }
    }

    fixedValueFvPatchField<scalar>::updateCoeffs();

    // TODO: perform averaging for cells sharing more than one boundary face
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    fWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
