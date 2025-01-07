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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace clouds
{
    defineTypeNameAndDebug(parcel, 0);
    addToRunTimeSelectionTable(cloud, parcel, polyMesh);
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

void Foam::clouds::parcel::initialise(const bool predict)
{
    cloud::initialise(predict);
    coupled::initialise(predict);
}


void Foam::clouds::parcel::partition()
{
    cloud::partition();
    coupled::partition();
}


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

        if (context != contextType::functionObject)
        {
            result = Lagrangianm::initDdt(dimVolume, rhoc(subMesh)) || result;
        }
    }

    result = Lagrangianm::initDdt(dimVolume, U, dUdt) || result;

    if (context != contextType::functionObject)
    {
        result = Lagrangianm::initDdt(dimVolume, Uc(subMesh)) || result;
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

    // Solve the number equation if a model provides a number source
    if (LagrangianModels().addsSupToField(number))
    {
        LagrangianEqn<scalar> numberEqn
        (
            Lagrangianm::Ddt(deltaT, number)
         ==
            number*LagrangianModels().source(deltaT)
        );

        numberEqn.solve(final);
    }

    // Solve the mass equation if a model provides a mass source
    if (LagrangianModels().addsSupToField(m))
    {
        LagrangianEqn<scalar> mEqn
        (
            Lagrangianm::Ddt(deltaT, number, m)
         ==
            number*LagrangianModels().source(deltaT, m)
        );

        mEqn.solve(final);

        // Correct the diameter, assuming the density remains constant
        spherical::correct(toSubField(m/rho));

        if (context != contextType::functionObject && final)
        {
            LagrangianEqn<scalar> mcEqn
            (
                Lagrangianm::noDdt(deltaT, dimVolume, rhoc(subMesh))
             ==
                LagrangianModels().sourceProxy(deltaT, m, rhoc(subMesh))
            );

            carrierEqn(rhoc) += number*mcEqn;
        }
    }

    // Solve the velocity equation
    {
        LagrangianEqn<vector> UEqn
        (
            Lagrangianm::Ddt(deltaT, m, U)
         ==
            LagrangianModels().source(deltaT, m, U)
        );

        UEqn.solve(final);

        if (context != contextType::functionObject && final)
        {
            LagrangianEqn<vector> UcEqn
            (
                Lagrangianm::noDdt(deltaT, dimMass, Uc(subMesh))
             ==
                LagrangianModels().sourceProxy(deltaT, m, U, Uc(subMesh))
            );

            carrierEqn(Uc) += number*UcEqn;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::clouds::parcel::parcel
(
    const polyMesh& mesh,
    const word& name,
    const contextType context,
    const IOobject::readOption readOption,
    const IOobject::writeOption writeOption
)
:
    cloud(mesh, name, context, readOption, writeOption),
    grouped(static_cast<const cloud&>(*this)),
    spherical(*this, *this),
    massive(*this, *this),
    coupledToFluid(static_cast<const cloud&>(*this)),
    sphericalCoupled(*this, *this, *this)
{
    reCalculateModified();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::clouds::parcel::~parcel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::clouds::parcel::solve()
{
    // Pre-solve operations ...
    coupled::clearCarrierEqns();

    cloud::solve();

    // Post-solve operations ...
}


// ************************************************************************* //
