/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "pressureWork.H"
#include "addToRunTimeSelectionTable.H"
#include "coupledToThermalFluid.H"
#include "fluidLagrangianThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace Lagrangian
{
    defineTypeNameAndDebug(pressureWork, 0);
    addToRunTimeSelectionTable(LagrangianModel, pressureWork, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Lagrangian::pressureWork::pressureWork
(
    const word& name,
    const LagrangianMesh& mesh,
    const dictionary& modelDict,
    const dictionary& stateDict
)
:
    LagrangianModel(name, mesh),
    cloudLagrangianModel(static_cast<const LagrangianModel&>(*this)),
    thermalCloud_(cloud<clouds::thermal>()),
    pPrevPtr_(nullptr)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::wordList Foam::Lagrangian::pressureWork::addSupFields() const
{
    return wordList({thermalCloud_.e.name()});
}


bool Foam::Lagrangian::pressureWork::addsSupToField
(
    const word& fieldName,
    const word& eqnFieldName
) const
{
    return
        fieldName == thermalCloud_.e.name()
     && (
            eqnFieldName == thermalCloud_.e.name()
         || eqnFieldName == cloud<clouds::coupledToThermalFluid>().hec.name()
        );
}


void Foam::Lagrangian::pressureWork::preAddSup
(
    const LagrangianSubScalarField& deltaT,
    const bool final
)
{
    const LagrangianSubMesh& subMesh = deltaT.mesh();

    const fluidLagrangianThermo& fluidThermo =
        thermalCloud_.thermo<fluidLagrangianThermo>(*this);
    const LagrangianSubScalarSubField p(subMesh.sub(fluidThermo.p()));

    pPrevPtr_.set(new LagrangianSubScalarField(p));
}


void Foam::Lagrangian::pressureWork::addSup
(
    const LagrangianSubScalarField& deltaT,
    const LagrangianSubScalarSubField& vOrM,
    const LagrangianSubScalarSubField& e,
    LagrangianEqn<scalar>& eqn
) const
{
    const LagrangianSubMesh& subMesh = deltaT.mesh();

    const LagrangianSubScalarField& m = thermalCloud_.m(subMesh);
    const LagrangianSubScalarSubField& rho = thermalCloud_.rho(subMesh);

    const fluidLagrangianThermo& fluidThermo =
        thermalCloud_.thermo<fluidLagrangianThermo>(*this);
    const LagrangianSubScalarSubField p(subMesh.sub(fluidThermo.p()));

    // Evaluate the time derivative of the specific volume
    tmp<LagrangianEqn<scalar>> dOneByRhoDt =
        - Lagrangianm::ddt0(rho)/(rho*rho.oldTime());

    // If this is the particle energy equation, then construct corrections to
    // the above time derivative. If not then just apply it directly.
    if (eqn.isPsi(e))
    {
        tmp<LagrangianEqn<scalar>> dOneByRhoDtCorr =
            tmp<LagrangianEqn<scalar>>
            (
                new LagrangianEqn<scalar>(subMesh)
            );

        // Correct for changes in pressure
        const LagrangianSubScalarSubField psi(subMesh.sub(fluidThermo.psi()));
        dOneByRhoDtCorr.ref().deltaTSu -=
            psi*(p - pPrevPtr_())/(rho*rho.oldTime());

        // Correct for changes in energy
        const LagrangianSubScalarSubField Cv(subMesh.sub(fluidThermo.Cv()));
        tmp<LagrangianSubScalarField> alphavByRhoCv =
            fluidThermo.alphav(subMesh)/(rho*Cv);
        dOneByRhoDtCorr.ref().deltaTSp += alphavByRhoCv();
        dOneByRhoDtCorr.ref().deltaTSu -= alphavByRhoCv()*e;

        // Correct for non-ideal thermodynamic modelling
        tmp<LagrangianSubScalarField> gamma(fluidThermo.Cp(subMesh)/Cv);
        dOneByRhoDtCorr.ref() /= gamma - p*alphavByRhoCv;

        eqn -= m*p*(dOneByRhoDt + dOneByRhoDtCorr);
    }
    else
    {
        const LagrangianSubScalarField& pc =
            cloud<clouds::coupledToThermalFluid>().pc(subMesh);

        eqn -= m*pc*dOneByRhoDt;
    }
}


void Foam::Lagrangian::pressureWork::postAddSup
(
    const LagrangianSubScalarField& deltaT,
    const bool final
)
{
    pPrevPtr_.clear();
}


// ************************************************************************* //
