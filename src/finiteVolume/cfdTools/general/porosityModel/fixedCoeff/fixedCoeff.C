/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2024 OpenFOAM Foundation
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

#include "addToRunTimeSelectionTable.H"
#include "fixedCoeff.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace porosityModels
    {
        defineTypeNameAndDebug(fixedCoeff, 0);
        addToRunTimeSelectionTable(porosityModel, fixedCoeff, mesh);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::porosityModels::fixedCoeff::apply
(
    scalarField& Udiag,
    vectorField& Usource,
    const scalarField& V,
    const vectorField& U,
    const scalar rho
) const
{
    const labelList& cells = mesh_.cellZones()[zoneName_];

    forAll(cells, i)
    {
        const label celli = cells[i];
        const label j = fieldIndex(i);
        const tensor Cd = rho*(alpha_[j] + beta_[j]*mag(U[celli]));
        const scalar isoCd = tr(Cd);

        Udiag[celli] += V[celli]*isoCd;
        Usource[celli] -= V[celli]*((Cd - I*isoCd) & U[celli]);
    }
}


void Foam::porosityModels::fixedCoeff::apply
(
    tensorField& AU,
    const vectorField& U,
    const scalar rho
) const
{
    const labelList& cells = mesh_.cellZones()[zoneName_];

    forAll(cells, i)
    {
        const label celli = cells[i];
        const label j = fieldIndex(i);
        const tensor alpha = alpha_[j];
        const tensor beta = beta_[j];

        AU[celli] += rho*(alpha + beta*mag(U[celli]));
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porosityModels::fixedCoeff::fixedCoeff
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& coeffDict,
    const word& cellZoneName
)
:
    porosityModel(name, mesh, dict, coeffDict, cellZoneName),
    alphaXYZ_("alpha", dimless/dimTime, coeffDict),
    betaXYZ_("beta", dimless/dimLength, coeffDict),
    rhoRefFound_(coeffDict.found("rhoRef")),
    rhoRef_(coeffDict.lookupOrDefault<scalar>("rhoRef", 1.0))
{
    adjustNegativeResistance(alphaXYZ_);
    adjustNegativeResistance(betaXYZ_);

    calcTransformModelData();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::porosityModels::fixedCoeff::~fixedCoeff()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::porosityModels::fixedCoeff::calcTransformModelData()
{
    if (coordSys_.R().uniform())
    {
        alpha_.setSize(1);
        beta_.setSize(1);

        alpha_[0] = Zero;
        alpha_[0].xx() = alphaXYZ_.value().x();
        alpha_[0].yy() = alphaXYZ_.value().y();
        alpha_[0].zz() = alphaXYZ_.value().z();
        alpha_[0] = coordSys_.R().transform(Zero, alpha_[0]);

        beta_[0] = Zero;
        beta_[0].xx() = betaXYZ_.value().x();
        beta_[0].yy() = betaXYZ_.value().y();
        beta_[0].zz() = betaXYZ_.value().z();
        beta_[0] = coordSys_.R().transform(Zero, beta_[0]);
    }
    else
    {
        const labelList& cells = mesh_.cellZones()[zoneName_];

        alpha_.setSize(cells.size());
        beta_.setSize(cells.size());

        forAll(cells, i)
        {
            alpha_[i] = Zero;
            alpha_[i].xx() = alphaXYZ_.value().x();
            alpha_[i].yy() = alphaXYZ_.value().y();
            alpha_[i].zz() = alphaXYZ_.value().z();

            beta_[i] = Zero;
            beta_[i].xx() = betaXYZ_.value().x();
            beta_[i].yy() = betaXYZ_.value().y();
            beta_[i].zz() = betaXYZ_.value().z();
        }

        const coordinateRotation& R = coordSys_.R
        (
            UIndirectList<vector>(mesh_.C(), cells)()
        );

        alpha_ = R.transform(alpha_);
        beta_ = R.transform(beta_);
    }
}


void Foam::porosityModels::fixedCoeff::calcForce
(
    const volVectorField& U,
    const volScalarField& rho,
    const volScalarField& mu,
    vectorField& force
) const
{
    scalarField Udiag(U.size(), 0.0);
    vectorField Usource(U.size(), Zero);
    const scalarField& V = mesh_.V();

    if (!rhoRefFound_)
    {
        FatalErrorInFunction
            << "rhoRef not specified" << exit(FatalError);
    }

    apply(Udiag, Usource, V, U, rhoRef_);

    force = Udiag*U - Usource;
}


void Foam::porosityModels::fixedCoeff::correct
(
    fvVectorMatrix& UEqn
) const
{
    const vectorField& U = UEqn.psi();
    const scalarField& V = mesh_.V();
    scalarField& Udiag = UEqn.diag();
    vectorField& Usource = UEqn.source();

    if (UEqn.dimensions() == dimForce && !rhoRefFound_)
    {
        FatalErrorInFunction
            << "rhoRef not specified" << exit(FatalError);
    }

    apply(Udiag, Usource, V, U, rhoRef_);
}


void Foam::porosityModels::fixedCoeff::correct
(
    const fvVectorMatrix& UEqn,
    volTensorField& AU
) const
{
    const vectorField& U = UEqn.psi();

    if (UEqn.dimensions() == dimForce && !rhoRefFound_)
    {
        FatalErrorInFunction
            << "rhoRef not specified" << exit(FatalError);
    }

    apply(AU, U, rhoRef_);
}


// ************************************************************************* //
