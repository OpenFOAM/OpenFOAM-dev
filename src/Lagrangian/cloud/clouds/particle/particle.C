/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025-2026 OpenFOAM Foundation
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

#include "particle.H"
#include "cloud_fvModel.H"
#include "cloud_functionObject.H"
#include "LagrangiancDdt.H"
#include "LagrangianmDdt.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace clouds
{
    defineTypeNameAndDebug(particle, 0);
    addToRunTimeSelectionTable(cloud, particle, LagrangianMesh);
}
namespace fv
{
    makeCloudFvModel(particle);
}
namespace functionObjects
{
    makeCloudFunctionObject(particle);
}
}


// * * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * //

Foam::tmp<Foam::LagrangianSubVectorField> Foam::clouds::particle::dUdt
(
    const LagrangianSubMesh& subMesh
) const
{
    const LagrangianSubScalarSubField& m = this->m.ref(subMesh);
    const LagrangianSubVectorSubField& U = this->U.ref(subMesh);

    return
        LagrangianModels().addsSupToField(m)
      ? (Lagrangianc::Ddt(m, U) - Lagrangianc::Ddt(m)*U)/m
      : Lagrangianc::Ddt(U);
}


bool Foam::clouds::particle::reCalculateModified()
{
    const bool dUdt = tracking == trackingType::parabolic;

    const LagrangianSubMesh subMesh = this->mesh().subNone();

    LagrangianSubScalarSubField& m = this->m.ref(subMesh);
    LagrangianSubScalarSubField& e = this->e.ref(subMesh);
    LagrangianSubVectorSubField& U = this->U.ref(subMesh);

    bool result = false;

    if (LagrangianModels().addsSupToField(m))
    {
        result = Lagrangianm::initDdt(dimless, m, dUdt) || result;

        if (context == cloud::contextType::fvModel)
        {
            result = initPsicDdt(m, rhoc) || result;
            if (hasPhase())
            {
                result = initPsicDdt(m, rhocPhase) || result;
            }
        }
    }

    {
        result = Lagrangianm::initDdt(dimMass, e, dUdt) || result;

        if (context == cloud::contextType::fvModel)
        {
            result = initPsicDdt(m, hec) || result;
            if (hasPhase() && &hecPhase != &hec)
            {
                result = initPsicDdt(m, hecPhase) || result;
            }
        }
    }

    {
        result = Lagrangianm::initDdt(dimMass, U, dUdt) || result;

        if (context == cloud::contextType::fvModel)
        {
            result = initPsicDdt(m, Uc) || result;
            if (hasPhase() && &UcPhase != &Uc)
            {
                result = initPsicDdt(m, UcPhase) || result;
            }
        }
    }

    return result;
}


void Foam::clouds::particle::calculate
(
    const LagrangianSubScalarField& deltaT,
    const bool final
)
{
    const LagrangianSubMesh& subMesh = deltaT.mesh();

    LagrangianSubScalarSubField& m = this->m.ref(subMesh);
    const LagrangianSubScalarSubField& rho = this->rho(subMesh);
    LagrangianSubScalarSubField& e = this->e.ref(subMesh);
    LagrangianSubVectorSubField& U = this->U.ref(subMesh);

    // Update the pressure
    thermo().correctPressure(subMesh);

    // Solve the mass equation if a model provides a mass source
    if (LagrangianModels().addsSupToField(m))
    {
        LagrangianEqn<scalar> mEqn
        (
            Lagrangianm::Ddt(deltaT, m)
         ==
            LagrangianModels().source(deltaT, m)
        );

        mEqn.solve(final);

        // Correct the diameter for the change in mass, assuming the density
        // remains constant
        spherical::correct(toSubField(m/rho));

        // Calculate mass exchanges with the carrier
        if (context == cloud::contextType::fvModel && final)
        {
            carrierEqn(rhoc) += psicEqn(deltaT, m, rhoc);
            if (hasPhase())
            {
                carrierEqn(rhocPhase) += psicEqn(deltaT, m, rhocPhase);
            }
        }
    }

    // Solve the energy equation
    {
        LagrangianEqn<scalar> eEqn
        (
            Lagrangianm::Ddt(deltaT, m, e)
         ==
            LagrangianModels().source(deltaT, m, e)
        );

        eEqn.solve(final);

        // Update the thermodynamic model
        thermo().correct(subMesh);

        // Correct the diameter for changes in density
        spherical::correct(toSubField(m/rho));

        // Calculate energy exchanges with the carrier
        if (context == cloud::contextType::fvModel && final)
        {
            carrierEqn(hec) += psicEqn(deltaT, m, e, hec);
            if (hasPhase() && &hecPhase != &hec)
            {
                carrierEqn(hecPhase) += psicEqn(deltaT, m, e, hecPhase);
            }
        }
    }

    // Solve the momentum equation
    {
        LagrangianEqn<vector> UEqn
        (
            Lagrangianm::Ddt(deltaT, m, U)
         ==
            LagrangianModels().source(deltaT, m, U)
        );

        UEqn.solve(final);

        // Calculate momentum exchanges with the carrier
        if (context == cloud::contextType::fvModel && final)
        {
            carrierEqn(Uc) += psicEqn(deltaT, m, U, Uc);
            if (hasPhase() && &UcPhase != &Uc)
            {
                carrierEqn(UcPhase) += psicEqn(deltaT, m, U, UcPhase);
            }
        }
    }
}


void Foam::clouds::particle::partition()
{
    cloud::partition();
    carried::clearCarrierFields();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::clouds::particle::particle
(
    LagrangianMesh& mesh,
    const contextType context,
    const dictionary& dict
)
:
    cloud(mesh, context),
    spherical(static_cast<const cloud&>(*this)),
    coupledToThermalFluid(static_cast<const cloud&>(*this), dict),
    thermal(*this, *this, *this),
    sphericalCoupled(*this, *this, *this),
    massiveCoupledToFluid(*this, *this, *this)
{
    thermo().initialise();

    reCalculateModified();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::clouds::particle::~particle()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::clouds::particle::solve(const bool initial, const bool final)
{
    // Pre-solve operations ...
    carried::resetCarrierFields(initial);
    coupled::clearCarrierEqns();
    coupledToThermalFluid::updateCarrier();

    // Solve
    cloud::solve(initial, final);

    // Post-solve operations ...
}


// ************************************************************************* //
