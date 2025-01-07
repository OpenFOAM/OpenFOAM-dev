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

#include "kinematicParcel.H"
#include "cloud_fvModel.H"
#include "cloud_functionObject.H"
#include "LagrangiancDdt.H"
#include "LagrangianmDdt.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace clouds
{
    defineTypeNameAndDebug(kinematicParcel, 0);
    addToRunTimeSelectionTable(cloud, kinematicParcel, polyMesh);
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

void Foam::clouds::kinematicParcel::initialise(const bool predict)
{
    cloud::initialise(predict);
    coupled::initialise(predict);
}


void Foam::clouds::kinematicParcel::partition()
{
    cloud::partition();
    coupled::partition();
}


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

        /*
        if (context != contextType::functionObject)
        {
            result = Lagrangianm::initDdt(dimVolume, onec()) || result;
        }
        */
    }

    result = Lagrangianm::initDdt(dimVolume, U, dUdt) || result;

    if (context != contextType::functionObject)
    {
        result = Lagrangianm::initDdt(dimVolume, Uc(subMesh)) || result;
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

    // Solve the number equation if a model provides a number source
    if (LagrangianModels().addsSupToField(word::null))
    {
        LagrangianEqn<scalar> numberEqn
        (
            Lagrangianm::Ddt(deltaT, number)
         ==
            number*LagrangianModels().source(deltaT)
        );

        numberEqn.solve(final);
    }

    // Solve the volume equation if a model provides a volume source
    if (LagrangianModels().addsSupToField(v))
    {
        LagrangianEqn<scalar> vEqn
        (
            Lagrangianm::Ddt(deltaT, v)
         ==
            LagrangianModels().source(deltaT, v)
        );

        vEqn.solve(final);

        // Correct the diameter
        spherical::correct(v);

        /*
        if (context != contextType::functionObject && final)
        {
            LagrangianEqn<scalar> vcEqn
            (
                Lagrangianm::noDdt(deltaT, dimVolume, onec())
             ==
                LagrangianModels().sourceProxy(deltaT, v, onec())
            );

            carrierEqn(one()) += number*vcEqn;
        }
        */
    }

    // Solve the velocity equation
    {
        LagrangianEqn<vector> UEqn
        (
            Lagrangianm::Ddt(deltaT, v, U)
         ==
            LagrangianModels().source(deltaT, v, U)
        );

        UEqn.solve(final);

        if (context != contextType::functionObject && final)
        {
            LagrangianEqn<vector> UcEqn
            (
                Lagrangianm::noDdt(deltaT, dimVolume, Uc(subMesh))
             ==
                LagrangianModels().sourceProxy(deltaT, v, U, Uc(subMesh))
            );

            carrierEqn(Uc) += number*UcEqn;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::clouds::kinematicParcel::kinematicParcel
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
    coupledToIncompressibleFluid(static_cast<const cloud&>(*this)),
    sphericalCoupled(*this, *this, *this)
{
    reCalculateModified();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::clouds::kinematicParcel::~kinematicParcel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::clouds::kinematicParcel::solve()
{
    // Pre-solve operations ...
    coupled::clearCarrierEqns();
    coupledToIncompressibleFluid::updateNuc();

    cloud::solve();

    // Post-solve operations ...
}


// ************************************************************************* //
