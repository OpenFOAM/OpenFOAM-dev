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

#include "kinematicParticle.H"
#include "cloud_fvModel.H"
#include "cloud_functionObject.H"
#include "LagrangiancDdt.H"
#include "LagrangianmDdt.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace clouds
{
    defineTypeNameAndDebug(kinematicParticle, 0);
    addToRunTimeSelectionTable(cloud, kinematicParticle, polyMesh);
}
namespace fv
{
    makeCloudFvModel(kinematicParticle);
}
namespace functionObjects
{
    makeCloudFunctionObject(kinematicParticle);
}
}


// * * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * //

void Foam::clouds::kinematicParticle::initialise(const bool predict)
{
    cloud::initialise(predict);
    coupled::initialise(predict);
}


void Foam::clouds::kinematicParticle::partition()
{
    cloud::partition();
    coupled::partition();
}


Foam::tmp<Foam::LagrangianSubVectorField>
Foam::clouds::kinematicParticle::dUdt
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


bool Foam::clouds::kinematicParticle::reCalculateModified()
{
    const bool dUdt = tracking == trackingType::parabolic;

    const LagrangianSubMesh subMesh = this->mesh().subNone();

    LagrangianSubScalarSubField& v = this->v.ref(subMesh);
    LagrangianSubVectorSubField& U = this->U.ref(subMesh);

    bool result = false;

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


void Foam::clouds::kinematicParticle::calculate
(
    const LagrangianSubScalarField& deltaT,
    const bool final
)
{
    const LagrangianSubMesh& subMesh = deltaT.mesh();

    LagrangianSubScalarSubField& v = this->v.ref(subMesh);
    LagrangianSubVectorSubField& U = this->U.ref(subMesh);

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

            carrierEqn(one()) += vcEqn;
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

            carrierEqn(Uc) += UcEqn;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::clouds::kinematicParticle::kinematicParticle
(
    const polyMesh& mesh,
    const word& name,
    const contextType context,
    const IOobject::readOption readOption,
    const IOobject::writeOption writeOption
)
:
    cloud(mesh, name, context, readOption, writeOption),
    spherical(static_cast<const cloud&>(*this)),
    coupledToIncompressibleFluid(static_cast<const cloud&>(*this)),
    sphericalCoupled(*this, *this, *this)
{
    reCalculateModified();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::clouds::kinematicParticle::~kinematicParticle()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::clouds::kinematicParticle::solve()
{
    // Pre-solve operations ...
    coupled::clearCarrierEqns();
    coupledToIncompressibleFluid::updateNuc();

    cloud::solve();

    // Post-solve operations ...
}


// ************************************************************************* //
