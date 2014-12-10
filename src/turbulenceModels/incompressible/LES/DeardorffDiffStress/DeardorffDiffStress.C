/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
namespace incompressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(DeardorffDiffStress, 0);
addToRunTimeSelectionTable(LESModel, DeardorffDiffStress, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void DeardorffDiffStress::updateSubGridScaleFields(const volScalarField& K)
{
    nuSgs_ = ck_*sqrt(K)*delta();
    nuSgs_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DeardorffDiffStress::DeardorffDiffStress
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    LESModel(modelName, U, phi, transport, turbulenceModelName),
    GenSGSStress(U, phi, transport),

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

    const volSymmTensorField D(symm(gradU));

    const volSymmTensorField P(-twoSymm(B_ & gradU));

    volScalarField K(0.5*tr(B_));
    volScalarField Epsilon(2*nuEff()*magSqr(D));

    tmp<fvSymmTensorMatrix> BEqn
    (
        fvm::ddt(B_)
      + fvm::div(phi(), B_)
      - fvm::laplacian(DBEff(), B_)
      + fvm::Sp(cm_*sqrt(K)/delta(), B_)
     ==
        P
      + 0.8*K*D
      - (2*ce_ - 0.667*cm_)*I*Epsilon
    );

    BEqn().relax();
    BEqn().solve();

    // Bounding the component kinetic energies

    forAll(B_, celli)
    {
        B_[celli].component(symmTensor::XX) =
            max(B_[celli].component(symmTensor::XX), kMin_.value());
        B_[celli].component(symmTensor::YY) =
            max(B_[celli].component(symmTensor::YY), kMin_.value());
        B_[celli].component(symmTensor::ZZ) =
            max(B_[celli].component(symmTensor::ZZ), kMin_.value());
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
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
