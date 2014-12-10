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

#include "DeardorffDiffStress.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(DeardorffDiffStress, 0);
addToRunTimeSelectionTable(LESModel, DeardorffDiffStress, dictionary);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void DeardorffDiffStress::updateSubGridScaleFields(const volScalarField& K)
{
    muSgs_ = ck_*rho()*sqrt(K)*delta();
    muSgs_.correctBoundaryConditions();

    alphaSgs_ = muSgs_/Prt_;
    alphaSgs_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DeardorffDiffStress::DeardorffDiffStress
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const fluidThermo& thermoPhysicalModel,
    const word& turbulenceModelName,
    const word& modelName
)
:
    LESModel(modelName, rho, U, phi, thermoPhysicalModel, turbulenceModelName),
    GenSGSStress(rho, U, phi, thermoPhysicalModel),

    ck_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ck",
            coeffDict_,
            0.094
        )
    ),
    cm_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cm",
            coeffDict_,
            4.13
        )
    )
{
    updateSubGridScaleFields(0.5*tr(B_));

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void DeardorffDiffStress::correct(const tmp<volTensorField>& tgradU)
{
    const volTensorField& gradU = tgradU();

    GenSGSStress::correct(gradU);

    volSymmTensorField D(symm(gradU));

    volSymmTensorField P(-rho()*twoSymm(B_ & gradU));

    volScalarField K(0.5*tr(B_));

    tmp<fvSymmTensorMatrix> BEqn
    (
        fvm::ddt(rho(), B_)
      + fvm::div(phi(), B_)
      - fvm::laplacian(DBEff(), B_)
      + fvm::Sp(cm_*rho()*sqrt(K)/delta(), B_)
     ==
        P
      + 0.8*rho()*K*D
      - (2*ce_ - 0.667*cm_)*I*rho()*epsilon()
    );

    BEqn().relax();
    BEqn().solve();

    // Bounding the component kinetic energies

    forAll(B_, celli)
    {
        B_[celli].component(symmTensor::XX) =
            max(B_[celli].component(symmTensor::XX), 1.0e-10);
        B_[celli].component(symmTensor::YY) =
            max(B_[celli].component(symmTensor::YY), 1.0e-10);
        B_[celli].component(symmTensor::ZZ) =
            max(B_[celli].component(symmTensor::ZZ), 1.0e-10);
    }

    K = 0.5*tr(B_);
    bound(K, kMin_);

    updateSubGridScaleFields(K);
}


bool DeardorffDiffStress::read()
{
    if (GenSGSStress::read())
    {
        ck_.readIfPresent(coeffDict());
        cm_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
