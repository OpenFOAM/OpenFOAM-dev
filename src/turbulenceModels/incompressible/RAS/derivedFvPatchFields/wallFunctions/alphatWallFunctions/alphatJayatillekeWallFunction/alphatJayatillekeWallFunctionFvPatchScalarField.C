/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
#include "incompressible/turbulenceModel/turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

scalar alphatJayatillekeWallFunctionFvPatchScalarField::tolerance_ = 0.01;
label alphatJayatillekeWallFunctionFvPatchScalarField::maxIters_ = 10;

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void alphatJayatillekeWallFunctionFvPatchScalarField::checkType()
{
    if (!isA<wallFvPatch>(patch()))
    {
        FatalErrorIn
        (
            "alphatJayatillekeWallFunctionFvPatchScalarField::checkType()"
        )   << "Invalid wall function specification" << nl
            << "    Patch type for patch " << patch().name()
            << " must be wall" << nl
            << "    Current patch type is " << patch().type() << nl << endl
            << abort(FatalError);
    }
}


scalar alphatJayatillekeWallFunctionFvPatchScalarField::Psmooth
(
    const scalar Prat
) const
{
    return 9.24*(pow(Prat, 0.75) - 1.0)*(1.0 + 0.28*exp(-0.007*Prat));
}


scalar alphatJayatillekeWallFunctionFvPatchScalarField::yPlusTherm
(
    const scalar P,
    const scalar Prat
) const
{
    scalar ypt = 11.0;

    for (int i=0; i<maxIters_; i++)
    {
        scalar f = ypt - (log(E_*ypt)/kappa_ + P)/Prat;
        scalar df = 1.0 - 1.0/(ypt*kappa_*Prat);
        scalar yptNew = ypt - f/df;

        if (yptNew < VSMALL)
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
    Prt_(0.85),
    Cmu_(0.09),
    kappa_(0.41),
    E_(9.8)
{
    checkType();
}


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
    Prt_(ptf.Prt_),
    Cmu_(ptf.Cmu_),
    kappa_(ptf.kappa_),
    E_(ptf.E_)
{
    checkType();
}


alphatJayatillekeWallFunctionFvPatchScalarField::
alphatJayatillekeWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    Prt_(readScalar(dict.lookup("Prt"))), // force read to avoid ambiguity
    Cmu_(dict.lookupOrDefault<scalar>("Cmu", 0.09)),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    E_(dict.lookupOrDefault<scalar>("E", 9.8))
{
    checkType();
}


alphatJayatillekeWallFunctionFvPatchScalarField::
alphatJayatillekeWallFunctionFvPatchScalarField
(
    const alphatJayatillekeWallFunctionFvPatchScalarField& wfpsf
)
:
    fixedValueFvPatchScalarField(wfpsf),
    Prt_(wfpsf.Prt_),
    Cmu_(wfpsf.Cmu_),
    kappa_(wfpsf.kappa_),
    E_(wfpsf.E_)
{
    checkType();
}


alphatJayatillekeWallFunctionFvPatchScalarField::
alphatJayatillekeWallFunctionFvPatchScalarField
(
    const alphatJayatillekeWallFunctionFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(wfpsf, iF),
    Prt_(wfpsf.Prt_),
    Cmu_(wfpsf.Cmu_),
    kappa_(wfpsf.kappa_),
    E_(wfpsf.E_)
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void alphatJayatillekeWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchi = patch().index();

    // Retrieve turbulence properties from model
    const turbulenceModel& turbModel =
        db().lookupObject<turbulenceModel>("turbulenceModel");
    const scalar Cmu25 = pow(Cmu_, 0.25);
    const scalarField& y = turbModel.y()[patchi];
    const tmp<volScalarField> tnu = turbModel.nu();
    const volScalarField& nu = tnu();
    const scalarField& nuw = nu.boundaryField()[patchi];
    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();

    const IOdictionary& transportProperties =
        db().lookupObject<IOdictionary>("transportProperties");

    // Molecular Prandtl number
    const scalar Pr
    (
        dimensionedScalar(transportProperties.lookup("Pr")).value()
    );

    // Populate boundary values
    scalarField& alphatw = *this;
    forAll(alphatw, faceI)
    {
        label faceCellI = patch().faceCells()[faceI];

        // y+
        scalar yPlus = Cmu25*sqrt(k[faceCellI])*y[faceI]/nuw[faceI];

        // Molecular-to-turbulent Prandtl number ratio
        scalar Prat = Pr/Prt_;

        // Thermal sublayer thickness
        scalar P = Psmooth(Prat);
        scalar yPlusTherm = this->yPlusTherm(P, Prat);

        // Update turbulent thermal conductivity
        if (yPlus > yPlusTherm)
        {
            scalar nu = nuw[faceI];
            scalar kt = nu*(yPlus/(Prt_*(log(E_*yPlus)/kappa_ + P)) - 1/Pr);
            alphatw[faceI] = max(0.0, kt);
        }
        else
        {
            alphatw[faceI] = 0.0;
        }
    }

    fixedValueFvPatchField<scalar>::updateCoeffs();
}


void alphatJayatillekeWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeKeyword("Prt") << Prt_ << token::END_STATEMENT << nl;
    os.writeKeyword("Cmu") << Cmu_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappa") << kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("E") << E_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    alphatJayatillekeWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
