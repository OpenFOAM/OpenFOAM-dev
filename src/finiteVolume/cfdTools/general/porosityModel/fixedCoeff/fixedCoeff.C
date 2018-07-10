/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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
    forAll(cellZoneIDs_, zoneI)
    {
        const tensorField& alphaZones = alpha_[zoneI];
        const tensorField& betaZones = beta_[zoneI];

        const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];

        forAll(cells, i)
        {
            const label celli = cells[i];
            const label j = fieldIndex(i);
            const tensor Cd = rho*(alphaZones[j] + betaZones[j]*mag(U[celli]));
            const scalar isoCd = tr(Cd);

            Udiag[celli] += V[celli]*isoCd;
            Usource[celli] -= V[celli]*((Cd - I*isoCd) & U[celli]);
        }
    }
}


void Foam::porosityModels::fixedCoeff::apply
(
    tensorField& AU,
    const vectorField& U,
    const scalar rho
) const
{

    forAll(cellZoneIDs_, zoneI)
    {
        const tensorField& alphaZones = alpha_[zoneI];
        const tensorField& betaZones = beta_[zoneI];

        const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];

        forAll(cells, i)
        {
            const label celli = cells[i];
            const label j = fieldIndex(i);
            const tensor alpha = alphaZones[j];
            const tensor beta = betaZones[j];

            AU[celli] += rho*(alpha + beta*mag(U[celli]));
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porosityModels::fixedCoeff::fixedCoeff
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& cellZoneName
)
:
    porosityModel(name, modelType, mesh, dict, cellZoneName),
    alphaXYZ_("alpha", dimless/dimTime, coeffs_),
    betaXYZ_("beta", dimless/dimLength, coeffs_),
    alpha_(cellZoneIDs_.size()),
    beta_(cellZoneIDs_.size())
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
        forAll(cellZoneIDs_, zoneI)
        {
            alpha_[zoneI].setSize(1);
            beta_[zoneI].setSize(1);

            alpha_[zoneI][0] = Zero;
            alpha_[zoneI][0].xx() = alphaXYZ_.value().x();
            alpha_[zoneI][0].yy() = alphaXYZ_.value().y();
            alpha_[zoneI][0].zz() = alphaXYZ_.value().z();
            alpha_[zoneI][0] = coordSys_.R().transformTensor(alpha_[zoneI][0]);

            beta_[zoneI][0] = Zero;
            beta_[zoneI][0].xx() = betaXYZ_.value().x();
            beta_[zoneI][0].yy() = betaXYZ_.value().y();
            beta_[zoneI][0].zz() = betaXYZ_.value().z();
            beta_[zoneI][0] = coordSys_.R().transformTensor(beta_[zoneI][0]);
        }
    }
    else
    {
        forAll(cellZoneIDs_, zoneI)
        {
            const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];

            alpha_[zoneI].setSize(cells.size());
            beta_[zoneI].setSize(cells.size());

            forAll(cells, i)
            {
                alpha_[zoneI][i] = Zero;
                alpha_[zoneI][i].xx() = alphaXYZ_.value().x();
                alpha_[zoneI][i].yy() = alphaXYZ_.value().y();
                alpha_[zoneI][i].zz() = alphaXYZ_.value().z();

                beta_[zoneI][i] = Zero;
                beta_[zoneI][i].xx() = betaXYZ_.value().x();
                beta_[zoneI][i].yy() = betaXYZ_.value().y();
                beta_[zoneI][i].zz() = betaXYZ_.value().z();
            }

            const coordinateRotation& R = coordSys_.R(mesh_, cells);

            alpha_[zoneI] = R.transformTensor(alpha_[zoneI], cells);
            beta_[zoneI] = R.transformTensor(beta_[zoneI], cells);
        }
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
    scalar rhoRef = readScalar(coeffs_.lookup("rhoRef"));

    apply(Udiag, Usource, V, U, rhoRef);

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

    scalar rho = 1.0;
    if (UEqn.dimensions() == dimForce)
    {
        coeffs_.lookup("rhoRef") >> rho;
    }

    apply(Udiag, Usource, V, U, rho);
}


void Foam::porosityModels::fixedCoeff::correct
(
    fvVectorMatrix& UEqn,
    const volScalarField&,
    const volScalarField&
) const
{
    const vectorField& U = UEqn.psi();
    const scalarField& V = mesh_.V();
    scalarField& Udiag = UEqn.diag();
    vectorField& Usource = UEqn.source();

    scalar rho = 1.0;
    if (UEqn.dimensions() == dimForce)
    {
        coeffs_.lookup("rhoRef") >> rho;
    }

    apply(Udiag, Usource, V, U, rho);
}


void Foam::porosityModels::fixedCoeff::correct
(
    const fvVectorMatrix& UEqn,
    volTensorField& AU
) const
{
    const vectorField& U = UEqn.psi();

    scalar rho = 1.0;
    if (UEqn.dimensions() == dimForce)
    {
        coeffs_.lookup("rhoRef") >> rho;
    }

    apply(AU, U, rho);
}


bool Foam::porosityModels::fixedCoeff::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);

    return true;
}


// ************************************************************************* //
