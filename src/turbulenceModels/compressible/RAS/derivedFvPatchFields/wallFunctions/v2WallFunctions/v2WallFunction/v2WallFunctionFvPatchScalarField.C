/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2014 OpenFOAM Foundation
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

#include "v2WallFunctionFvPatchScalarField.H"
#include "RASModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void v2WallFunctionFvPatchScalarField::checkType()
{
    if (!isA<wallFvPatch>(patch()))
    {
        FatalErrorIn("v2WallFunctionFvPatchScalarField::checkType()")
            << "Invalid wall function specification" << nl
            << "    Patch type for patch " << patch().name()
            << " must be wall" << nl
            << "    Current patch type is " << patch().type() << nl << endl
            << abort(FatalError);
    }
}


void v2WallFunctionFvPatchScalarField::writeLocalEntries(Ostream& os) const
{
    os.writeKeyword("Cmu") << Cmu_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappa") << kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("E") << E_ << token::END_STATEMENT << nl;
}


scalar v2WallFunctionFvPatchScalarField::yPlusLam
(
    const scalar kappa,
    const scalar E
)
{
    scalar ypl = 11.0;

    for (int i=0; i<10; i++)
    {
        ypl = log(max(E*ypl, 1))/kappa;
    }

    return ypl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

v2WallFunctionFvPatchScalarField::v2WallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    Cmu_(0.09),
    kappa_(0.41),
    E_(9.8),
    yPlusLam_(yPlusLam(kappa_, E_))
{
    checkType();
}


v2WallFunctionFvPatchScalarField::v2WallFunctionFvPatchScalarField
(
    const v2WallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    Cmu_(ptf.Cmu_),
    kappa_(ptf.kappa_),
    E_(ptf.E_),
    yPlusLam_(ptf.yPlusLam_)
{
    checkType();
}


v2WallFunctionFvPatchScalarField::v2WallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
    Cmu_(dict.lookupOrDefault<scalar>("Cmu", 0.09)),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    E_(dict.lookupOrDefault<scalar>("E", 9.8)),
    yPlusLam_(yPlusLam(kappa_, E_))
{
    checkType();
}


v2WallFunctionFvPatchScalarField::v2WallFunctionFvPatchScalarField
(
    const v2WallFunctionFvPatchScalarField& v2wfpsf
)
:
    fixedValueFvPatchField<scalar>(v2wfpsf),
    Cmu_(v2wfpsf.Cmu_),
    kappa_(v2wfpsf.kappa_),
    E_(v2wfpsf.E_),
    yPlusLam_(v2wfpsf.yPlusLam_)
{
    checkType();
}


v2WallFunctionFvPatchScalarField::v2WallFunctionFvPatchScalarField
(
    const v2WallFunctionFvPatchScalarField& v2wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(v2wfpsf, iF),
    Cmu_(v2wfpsf.Cmu_),
    kappa_(v2wfpsf.kappa_),
    E_(v2wfpsf.E_),
    yPlusLam_(v2wfpsf.yPlusLam_)
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void v2WallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchI = patch().index();

    const turbulenceModel& turbulence =
        db().lookupObject<turbulenceModel>("turbulenceModel");
    const scalarField& y = turbulence.y()[patchI];

    const tmp<volScalarField> tk = turbulence.k();
    const volScalarField& k = tk();

    const tmp<volScalarField> tmu = turbulence.mu();
    const scalarField& muw = tmu().boundaryField()[patchI];

    const scalarField& rhow = turbulence.rho().boundaryField()[patchI];

    const scalar Cmu25 = pow025(Cmu_);

    scalarField& v2 = *this;

    // Set v2 wall values
    forAll(v2, faceI)
    {
        label faceCellI = patch().faceCells()[faceI];

        scalar uTau = Cmu25*sqrt(k[faceCellI]);

        scalar yPlus = uTau*y[faceI]/(muw[faceI]/rhow[faceI]);

        if (yPlus > yPlusLam_)
        {
            scalar Cv2 = 0.193;
            scalar Bv2 = -0.94;
            v2[faceI] = Cv2/kappa_*log(yPlus) + Bv2;
        }
        else
        {
            scalar Cv2 = 0.193;
            v2[faceI] = Cv2*pow4(yPlus);
        }

        v2[faceI] *= sqr(uTau);
    }

    fixedValueFvPatchField<scalar>::updateCoeffs();

    // TODO: perform averaging for cells sharing more than one boundary face
}


void v2WallFunctionFvPatchScalarField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    fixedValueFvPatchField<scalar>::evaluate(commsType);
}


void v2WallFunctionFvPatchScalarField::write(Ostream& os) const
{
    writeLocalEntries(os);
    fixedValueFvPatchField<scalar>::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    v2WallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
