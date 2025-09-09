/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "parcel.H"
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
    defineTypeNameAndDebug(parcel, 0);
    addToRunTimeSelectionTable(cloud, parcel, LagrangianMesh);
}
namespace fv
{
    makeCloudFvModel(parcel);
}
namespace functionObjects
{
    makeCloudFunctionObject(parcel);
}
}


// * * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * //

Foam::tmp<Foam::LagrangianSubVectorField> Foam::clouds::parcel::dUdt
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


bool Foam::clouds::parcel::reCalculateModified()
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

        result = reCalculateModified(m, rhoc) || result;

        if (hasPhase())
        {
            result = reCalculateModified(m, rhocPhase) || result;
        }
    }

    {
        result = Lagrangianm::initDdt(dimMass, U, dUdt) || result;

        result = reCalculateModified(m, Uc) || result;

        if (hasPhase() && &UcPhase != &Uc)
        {
            result = reCalculateModified(m, UcPhase) || result;
        }
    }

    return result;
}


void Foam::clouds::parcel::calculate
(
    const LagrangianSubScalarField& deltaT,
    const bool final
)
{
    const LagrangianSubMesh& subMesh = deltaT.mesh();

    LagrangianSubScalarSubField& number = this->number.ref(subMesh);
    LagrangianSubScalarSubField& m = this->m.ref(subMesh);
    const LagrangianSubScalarSubField rho(subMesh.sub(this->rho));
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

        calculate(deltaT, final, m, rhoc, number);

        if (hasPhase())
        {
            calculate(deltaT, final, m, rhocPhase, number);
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

        calculate(deltaT, final, m, U, Uc, number);

        if (hasPhase() && &UcPhase != &Uc)
        {
            calculate(deltaT, final, m, U, UcPhase, number);
        }
    }
}


void Foam::clouds::parcel::partition()
{
    cloud::partition();
    carried::clearCarrierFields();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::clouds::parcel::parcel
(
    LagrangianMesh& mesh,
    const contextType context,
    const dictionary& dict
)
:
    cloud(mesh, context),
    grouped(static_cast<const cloud&>(*this)),
    spherical(*this, *this),
    massive(*this, *this),
    coupledToFluid(static_cast<const cloud&>(*this), dict),
    sphericalCoupled(*this, *this, *this),
    massiveCoupledToFluid(*this, *this, *this)
{
    reCalculateModified();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::clouds::parcel::~parcel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::clouds::parcel::solve(const bool initial, const bool final)
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
