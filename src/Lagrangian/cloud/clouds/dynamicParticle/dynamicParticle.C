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

#include "dynamicParticle.H"
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
    defineTypeNameAndDebug(dynamicParticle, 0);
    addToRunTimeSelectionTable(cloud, dynamicParticle, LagrangianMesh);
}
namespace fv
{
    makeCloudFvModel(dynamicParticle);
}
namespace functionObjects
{
    makeCloudFunctionObject(dynamicParticle);
}
}


// * * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * //

Foam::tmp<Foam::LagrangianSubVectorField> Foam::clouds::dynamicParticle::dUdt
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


bool Foam::clouds::dynamicParticle::reCalculateModified()
{
    const bool dUdt = tracking == trackingType::parabolic;

    const LagrangianSubMesh subMesh = this->mesh().subNone();

    LagrangianSubScalarSubField& m = this->m.ref(subMesh);
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


void Foam::clouds::dynamicParticle::calculate
(
    const LagrangianSubScalarField& deltaT,
    const bool final
)
{
    const LagrangianSubMesh& subMesh = deltaT.mesh();

    LagrangianSubScalarSubField& m = this->m.ref(subMesh);
    const LagrangianSubScalarSubField& rho = this->rho(subMesh);
    LagrangianSubVectorSubField& U = this->U.ref(subMesh);

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

        // Correct the diameter, assuming the density remains constant
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


void Foam::clouds::dynamicParticle::partition()
{
    cloud::partition();
    carried::clearCarrierFields();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::clouds::dynamicParticle::dynamicParticle
(
    LagrangianMesh& mesh,
    const contextType context,
    const dictionary& dict
)
:
    cloud(mesh, context),
    spherical(static_cast<const cloud&>(*this)),
    coupledToFluid(static_cast<const cloud&>(*this), dict),
    dense(*this, *this),
    sphericalCoupled(*this, *this, *this),
    massiveCoupledToFluid(*this, *this, *this)
{
    reCalculateModified();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::clouds::dynamicParticle::~dynamicParticle()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::clouds::dynamicParticle::solve(const bool initial, const bool final)
{
    // Pre-solve operations ...
    carried::resetCarrierFields(initial);
    coupled::clearCarrierEqns();
    coupledToFluid::updateCarrier();

    // Solve
    cloud::solve(initial, final);

    // Post-solve operations ...
}


// ************************************************************************* //
