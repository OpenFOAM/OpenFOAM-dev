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

#include "nutkFilmWallFunctionFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fluidThermoMomentumTransportModel.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFilmRegionModel.H"
#include "mappedWallPolyPatch.H"
#include "mapDistribute.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField> nutkFilmWallFunctionFvPatchScalarField::calcUTau
(
    const scalarField& magGradU
) const
{
    tmp<scalarField> tuTau(new scalarField(patch().size(), 0.0));
    scalarField& uTau = tuTau.ref();

    typedef regionModels::surfaceFilmModels::surfaceFilmRegionModel modelType;

    bool foundFilm =
        db().time().foundObject<modelType>("surfaceFilmProperties");

    if (!foundFilm)
    {
        // Do nothing on construction - film model doesn't exist yet
        return tuTau;
    }

    const label patchi = patch().index();

    // Retrieve phase change mass from surface film model
    const modelType& filmModel =
        db().time().lookupObject<modelType>("surfaceFilmProperties");

    const label filmPatchi = filmModel.regionPatchID(patchi);

    tmp<volScalarField> mDotFilm(filmModel.primaryMassTrans());
    scalarField mDotFilmp = mDotFilm().boundaryField()[filmPatchi];
    filmModel.toPrimary(filmPatchi, mDotFilmp);


    // Retrieve RAS turbulence model
    const momentumTransportModel& turbModel =
        db().lookupObject<momentumTransportModel>
        (
            IOobject::groupName
            (
                momentumTransportModel::typeName,
                internalField().group()
            )
        );

    const scalarField& y = turbModel.y()[patchi];
    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const scalar Cmu25 = pow(Cmu_, 0.25);

    forAll(uTau, facei)
    {
        label faceCelli = patch().faceCells()[facei];

        scalar ut = Cmu25*sqrt(k[faceCelli]);

        scalar yPlus = y[facei]*ut/nuw[facei];

        scalar mStar = mDotFilmp[facei]/(y[facei]*ut);

        scalar factor = 0.0;
        if (yPlus > yPlusCrit_)
        {
            scalar expTerm = exp(min(50.0, B_*mStar));
            scalar powTerm = pow(yPlus, mStar/kappa_);
            factor = mStar/(expTerm*powTerm - 1.0 + rootVSmall);
        }
        else
        {
            scalar expTerm = exp(min(50.0, mStar));
            factor = mStar/(expTerm*yPlus - 1.0 + rootVSmall);
        }

        uTau[facei] = sqrt(max(0, magGradU[facei]*ut*factor));
    }

    return tuTau;
}


tmp<scalarField> nutkFilmWallFunctionFvPatchScalarField::nut() const
{
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

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalarField magGradU(mag(Uw.snGrad()));
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    return max
    (
        scalar(0),
        sqr(calcUTau(magGradU))/(magGradU + rootVSmall) - nuw
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nutkFilmWallFunctionFvPatchScalarField::nutkFilmWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutkWallFunctionFvPatchScalarField(p, iF),
    B_(5.5),
    yPlusCrit_(11.05)
{}


nutkFilmWallFunctionFvPatchScalarField::nutkFilmWallFunctionFvPatchScalarField
(
    const nutkFilmWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutkWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    B_(5.5),
    yPlusCrit_(11.05)
{}


nutkFilmWallFunctionFvPatchScalarField::nutkFilmWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutkWallFunctionFvPatchScalarField(p, iF, dict),
    B_(dict.lookupOrDefault("B", 5.5)),
    yPlusCrit_(dict.lookupOrDefault("yPlusCrit", 11.05))
{}


nutkFilmWallFunctionFvPatchScalarField::nutkFilmWallFunctionFvPatchScalarField
(
    const nutkFilmWallFunctionFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutkWallFunctionFvPatchScalarField(wfpsf, iF),
    B_(wfpsf.B_),
    yPlusCrit_(wfpsf.yPlusCrit_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> nutkFilmWallFunctionFvPatchScalarField::yPlus() const
{
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

    const scalarField& y = turbModel.y()[patchi];
    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    return y*calcUTau(mag(Uw.snGrad()))/nuw;
}


void nutkFilmWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    writeEntry(os, "B", B_);
    writeEntry(os, "yPlusCrit", yPlusCrit_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, nutkFilmWallFunctionFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
