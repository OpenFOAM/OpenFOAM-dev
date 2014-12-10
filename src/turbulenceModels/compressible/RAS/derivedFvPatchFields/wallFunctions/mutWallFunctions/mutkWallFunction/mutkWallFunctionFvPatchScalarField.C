/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "mutkWallFunctionFvPatchScalarField.H"
#include "compressible/turbulenceModel/turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField> mutkWallFunctionFvPatchScalarField::calcMut() const
{
    const label patchi = patch().index();
    const turbulenceModel& turbModel =
        db().lookupObject<turbulenceModel>("turbulenceModel");
    const scalarField& y = turbModel.y()[patchi];
    const scalarField& rhow = turbModel.rho().boundaryField()[patchi];
    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();
    const scalarField& muw = turbModel.mu().boundaryField()[patchi];

    const scalar Cmu25 = pow025(Cmu_);

    tmp<scalarField> tmutw(new scalarField(patch().size(), 0.0));
    scalarField& mutw = tmutw();

    forAll(mutw, faceI)
    {
        label faceCellI = patch().faceCells()[faceI];

        scalar yPlus =
            Cmu25*y[faceI]*sqrt(k[faceCellI])/(muw[faceI]/rhow[faceI]);

        if (yPlus > yPlusLam_)
        {
            mutw[faceI] = muw[faceI]*(yPlus*kappa_/log(E_*yPlus) - 1);
        }
    }

    return tmutw;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

mutkWallFunctionFvPatchScalarField::mutkWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mutWallFunctionFvPatchScalarField(p, iF)
{}


mutkWallFunctionFvPatchScalarField::mutkWallFunctionFvPatchScalarField
(
    const mutkWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mutWallFunctionFvPatchScalarField(ptf, p, iF, mapper)
{}


mutkWallFunctionFvPatchScalarField::mutkWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mutWallFunctionFvPatchScalarField(p, iF, dict)
{}


mutkWallFunctionFvPatchScalarField::mutkWallFunctionFvPatchScalarField
(
    const mutkWallFunctionFvPatchScalarField& wfpsf
)
:
    mutWallFunctionFvPatchScalarField(wfpsf)
{}


mutkWallFunctionFvPatchScalarField::mutkWallFunctionFvPatchScalarField
(
    const mutkWallFunctionFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mutWallFunctionFvPatchScalarField(wfpsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> mutkWallFunctionFvPatchScalarField::yPlus() const
{
    const label patchi = patch().index();

    const turbulenceModel& turbModel =
        db().lookupObject<turbulenceModel>("turbulenceModel");
    const scalarField& y = turbModel.y()[patchi];

    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();
    const scalarField kwc(k.boundaryField()[patchi].patchInternalField());
    const scalarField& muw = turbModel.mu().boundaryField()[patchi];
    const scalarField& rhow = turbModel.rho().boundaryField()[patchi];

    return pow025(Cmu_)*y*sqrt(kwc)/(muw/rhow);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    mutkWallFunctionFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
