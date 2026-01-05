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

#include "dynamicParcel.H"
#include "cloud_fvModel.H"
#include "cloud_functionObject.H"
#include "LagrangiancDdt.H"
#include "LagrangianmDdt.H"
#include "oneOrTmp.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace clouds
{
    defineTypeNameAndDebug(dynamicParcel, 0);
    addToRunTimeSelectionTable(cloud, dynamicParcel, LagrangianMesh);
}
namespace fv
{
    makeCloudFvModel(dynamicParcel);
}
namespace functionObjects
{
    makeCloudFunctionObject(dynamicParcel);
}
}


// * * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * //

Foam::tmp<Foam::LagrangianSubVectorField> Foam::clouds::dynamicParcel::dUdt
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


bool Foam::clouds::dynamicParcel::reCalculateModified()
{
    const bool dUdt = tracking == trackingType::parabolic;

    const LagrangianSubMesh subMesh = this->mesh().subNone();

    LagrangianSubScalarSubField& number = this->number.ref(subMesh);
    LagrangianSubScalarSubField& m = this->m.ref(subMesh);
    LagrangianSubVectorSubField& U = this->U.ref(subMesh);

    bool result = false;

    if (LagrangianModels().addsSupToField(word::null))
    {
        result = Lagrangianm::initDdt(dimless, number) || result;
    }

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


void Foam::clouds::dynamicParcel::calculate
(
    const LagrangianSubScalarField& deltaT,
    const bool final
)
{
    const LagrangianSubMesh& subMesh = deltaT.mesh();

    LagrangianSubScalarSubField& number = this->number.ref(subMesh);
    LagrangianSubScalarSubField& m = this->m.ref(subMesh);
    const LagrangianSubScalarSubField& rho = this->rho(subMesh);
    LagrangianSubVectorSubField& U = this->U.ref(subMesh);

    // Evaluate the fractional source
    LagrangianEqn<scalar> oneEqn(LagrangianModels().source(deltaT));

    // Initialise a unity fractional change in number (i.e., no change)
    oneOrTmp<LagrangianSubScalarField> numberByNumber0;

    // Solve the number equation if a model provides a fractional source
    if (oneEqn.valid())
    {
        LagrangianEqn<scalar> numberEqn
        (
            Lagrangianm::Ddt(deltaT, number)
         ==
            oneEqn
        );

        numberEqn.solve(final);

        // Set the fractional change in number
        numberByNumber0 = number/number.oldTime();

        // Correct the fractional source
        oneEqn *= numberByNumber0();
    }

    // Solve the mass equation if a model provides a mass source
    if (oneEqn.valid() || LagrangianModels().addsSupToField(m))
    {
        LagrangianEqn<scalar> mEqn
        (
            Lagrangianm::Ddt(deltaT, m)
          + oneEqn
         ==
            numberByNumber0()*LagrangianModels().source(deltaT, m)
        );

        mEqn.solve(final);

        // Correct the diameter, assuming the density remains constant
        spherical::correct(toSubField(m/rho));

        // Calculate mass exchanges with the carrier
        if (context == cloud::contextType::fvModel && final)
        {
            carrierEqn(rhoc) += number*psicEqn(deltaT, m, rhoc);
            if (hasPhase())
            {
                carrierEqn(rhocPhase) += number*psicEqn(deltaT, m, rhocPhase);
            }
        }
    }

    // Solve the velocity equation
    {
        LagrangianEqn<vector> UEqn
        (
            Lagrangianm::Ddt(deltaT, m, U)
          + m*oneEqn
         ==
            numberByNumber0*LagrangianModels().source(deltaT, m, U)
        );

        UEqn.solve(final);

        // Calculate momentum exchanges with the carrier
        if (context == cloud::contextType::fvModel && final)
        {
            carrierEqn(Uc) += number*psicEqn(deltaT, m, U, Uc);
            if (hasPhase() && &UcPhase != &Uc)
            {
                carrierEqn(UcPhase) += number*psicEqn(deltaT, m, U, UcPhase);
            }
        }
    }
}


void Foam::clouds::dynamicParcel::partition()
{
    cloud::partition();
    carried::clearCarrierFields();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::clouds::dynamicParcel::dynamicParcel
(
    LagrangianMesh& mesh,
    const contextType context,
    const dictionary& dict
)
:
    cloud(mesh, context),
    grouped(static_cast<const cloud&>(*this)),
    spherical(*this, *this),
    coupledToFluid(static_cast<const cloud&>(*this), dict),
    dense(*this, *this),
    sphericalCoupled(*this, *this, *this),
    massiveCoupledToFluid(*this, *this, *this)
{
    reCalculateModified();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::clouds::dynamicParcel::~dynamicParcel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::clouds::dynamicParcel::solve(const bool initial, const bool final)
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
