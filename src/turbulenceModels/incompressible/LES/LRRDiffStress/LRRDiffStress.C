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

#include "LRRDiffStress.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(LRRDiffStress, 0);
addToRunTimeSelectionTable(LESModel, LRRDiffStress, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void LRRDiffStress::updateSubGridScaleFields(const volScalarField& K)
{
    nuSgs_ = ck_*sqrt(K)*delta();
    nuSgs_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

LRRDiffStress::LRRDiffStress
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
            0.09
        )
    ),
    c1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "c1",
            coeffDict_,
            1.8
        )
    ),
    c2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "c2",
            coeffDict_,
            0.6
        )
    )
{
    volScalarField K(0.5*tr(B_));
    bound(K, kMin_);

    updateSubGridScaleFields(K);

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void LRRDiffStress::correct(const tmp<volTensorField>& tgradU)
{
    const volTensorField& gradU = tgradU();

    GenSGSStress::correct(gradU);

    const volSymmTensorField D(symm(gradU));

    const volSymmTensorField P(-twoSymm(B_ & gradU));

    volScalarField K(0.5*tr(B_));
    const volScalarField Epsilon(2*nuEff()*magSqr(D));

    tmp<fvSymmTensorMatrix> BEqn
    (
        fvm::ddt(B_)
      + fvm::div(phi(), B_)
      - fvm::laplacian(DBEff(), B_)
      + fvm::Sp(c1_*Epsilon/K, B_)
     ==
        P
      - (0.667*(1.0 - c1_)*I)*Epsilon
      - c2_*(P - 0.333*I*tr(P))
      - (0.667 - 2*c1_)*I*pow(K, 1.5)/delta()
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


bool LRRDiffStress::read()
{
    if (GenSGSStress::read())
    {
        ck_.readIfPresent(coeffDict());
        c1_.readIfPresent(coeffDict());
        c2_.readIfPresent(coeffDict());

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
