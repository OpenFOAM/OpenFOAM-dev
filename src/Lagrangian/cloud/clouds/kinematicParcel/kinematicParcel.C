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

#include "kinematicParcel.H"
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
    defineTypeNameAndDebug(kinematicParcel, 0);
    addToRunTimeSelectionTable(cloud, kinematicParcel, LagrangianMesh);
}
namespace fv
{
    makeCloudFvModel(kinematicParcel);
}
namespace functionObjects
{
    makeCloudFunctionObject(kinematicParcel);
}
}


// * * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * //

Foam::tmp<Foam::LagrangianSubVectorField>
Foam::clouds::kinematicParcel::dUdt
(
    const LagrangianSubMesh& subMesh
) const
{
    const LagrangianSubScalarSubField& v = this->v.ref(subMesh);
    const LagrangianSubVectorSubField& U = this->U.ref(subMesh);

    return
        LagrangianModels().addsSupToField(v)
      ? (Lagrangianc::Ddt(v, U) - Lagrangianc::Ddt(v)*U)/v
      : Lagrangianc::Ddt(U);
}


bool Foam::clouds::kinematicParcel::reCalculateModified()
{
    const bool dUdt = tracking == trackingType::parabolic;

    const LagrangianSubMesh subMesh = this->mesh().subNone();

    LagrangianSubScalarSubField& number = this->number.ref(subMesh);
    LagrangianSubScalarSubField& v = this->v.ref(subMesh);
    LagrangianSubVectorSubField& U = this->U.ref(subMesh);

    bool result = false;

    if (LagrangianModels().addsSupToField(word::null))
    {
        result = Lagrangianm::initDdt(dimless, number) || result;
    }

    if (LagrangianModels().addsSupToField(v))
    {
        result = Lagrangianm::initDdt(dimless, v, dUdt) || result;

        if (context == cloud::contextType::fvModel)
        {
            result = initPsicDdt(v, onec) || result;
            if (hasPhase())
            {
                result = initPsicDdt(v, onecPhase) || result;
            }
        }
    }

    {
        result = Lagrangianm::initDdt(dimVolume, U, dUdt) || result;

        if (context == cloud::contextType::fvModel)
        {
            result = initPsicDdt(v, Uc) || result;
            if (hasPhase() && &UcPhase != &Uc)
            {
                result = initPsicDdt(v, UcPhase) || result;
            }
        }
    }

    return result;
}


void Foam::clouds::kinematicParcel::calculate
(
    const LagrangianSubScalarField& deltaT,
    const bool final
)
{
    const LagrangianSubMesh& subMesh = deltaT.mesh();

    LagrangianSubScalarSubField& number = this->number.ref(subMesh);
    LagrangianSubScalarSubField& v = this->v.ref(subMesh);
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

    // Solve the volume equation if a model provides a volume source
    if (oneEqn.valid() || LagrangianModels().addsSupToField(v))
    {
        LagrangianEqn<scalar> vEqn
        (
            Lagrangianm::Ddt(deltaT, v)
          + oneEqn
         ==
            numberByNumber0()*LagrangianModels().source(deltaT, v)
        );

        vEqn.solve(final);

        // Correct the diameter
        spherical::correct(v);

        // Calculate volume exchanges with the carrier
        if (context == cloud::contextType::fvModel && final)
        {
            carrierEqn(onec) += number*psicEqn(deltaT, v, onec);
            if (hasPhase())
            {
                carrierEqn(onecPhase) += number*psicEqn(deltaT, v, onecPhase);
            }
        }
    }

    // Solve the velocity equation
    {
        LagrangianEqn<vector> UEqn
        (
            Lagrangianm::Ddt(deltaT, v, U)
          + v*oneEqn
         ==
            numberByNumber0*LagrangianModels().source(deltaT, v, U)
        );

        UEqn.solve(final);

        // Calculate momentum exchanges with the carrier
        if (context == cloud::contextType::fvModel && final)
        {
            carrierEqn(Uc) += number*psicEqn(deltaT, v, U, Uc);
            if (hasPhase() && &UcPhase != &Uc)
            {
                carrierEqn(UcPhase) += number*psicEqn(deltaT, v, U, UcPhase);
            }
        }
    }
}


void Foam::clouds::kinematicParcel::partition()
{
    cloud::partition();
    carried::clearCarrierFields();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::clouds::kinematicParcel::kinematicParcel
(
    LagrangianMesh& mesh,
    const contextType context,
    const dictionary& dict
)
:
    cloud(mesh, context),
    grouped(static_cast<const cloud&>(*this)),
    spherical(*this, *this),
    coupledToConstantDensityFluid(static_cast<const cloud&>(*this), dict),
    sphericalCoupled(*this, *this, *this)
{
    reCalculateModified();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::clouds::kinematicParcel::~kinematicParcel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::clouds::kinematicParcel::solve
(
    const bool initial,
    const bool final
)
{
    // Pre-solve operations ...
    carried::resetCarrierFields(initial);
    coupled::clearCarrierEqns();
    coupledToConstantDensityFluid::updateCarrier();

    // Solve
    cloud::solve(initial, final);

    // Post-solve operations ...
}


// ************************************************************************* //
