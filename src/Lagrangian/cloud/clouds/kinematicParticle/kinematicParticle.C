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
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace clouds
{
    defineTypeNameAndDebug(kinematicParticle, 0);
    addToRunTimeSelectionTable(cloud, kinematicParticle, LagrangianMesh);
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

        result = reCalculateModified(v, onec) || result;

        if (hasPhase())
        {
            result = reCalculateModified(v, onecPhase) || result;
        }
    }

    {
        result = Lagrangianm::initDdt(dimVolume, U, dUdt) || result;

        result = reCalculateModified(v, Uc) || result;

        if (hasPhase() && &UcPhase != &Uc)
        {
            result = reCalculateModified(v, UcPhase) || result;
        }
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

        calculate(deltaT, final, v, onec);

        if (hasPhase())
        {
            calculate(deltaT, final, v, onecPhase);
        }
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

        calculate(deltaT, final, v, U, Uc);

        if (hasPhase() && &UcPhase != &Uc)
        {
            calculate(deltaT, final, v, U, UcPhase);
        }
    }
}


void Foam::clouds::kinematicParticle::partition()
{
    cloud::partition();
    carried::clearCarrierFields();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::clouds::kinematicParticle::kinematicParticle
(
    LagrangianMesh& mesh,
    const contextType context,
    const dictionary& dict
)
:
    cloud(mesh, context),
    spherical(static_cast<const cloud&>(*this)),
    coupledToConstantDensityFluid(static_cast<const cloud&>(*this), dict),
    sphericalCoupled(*this, *this, *this)
{
    reCalculateModified();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::clouds::kinematicParticle::~kinematicParticle()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::clouds::kinematicParticle::solve
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
