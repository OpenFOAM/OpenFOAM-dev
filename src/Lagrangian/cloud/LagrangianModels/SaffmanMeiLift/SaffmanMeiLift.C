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

#include "SaffmanMeiLift.H"
#include "addToRunTimeSelectionTable.H"
#include "coupledToConstantDensityFluid.H"
#include "coupledToFluid.H"
#include "sphericalCoupled.H"
#include "LagrangianmSp.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace Lagrangian
{
    defineTypeNameAndDebug(SaffmanMeiLift, 0);
    addToRunTimeSelectionTable
    (
        LagrangianModel,
        SaffmanMeiLift,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::LagrangianSubTensorField>
Foam::Lagrangian::SaffmanMeiLift::calcL
(
    const LagrangianModelRef& model,
    const LagrangianSubMesh& subMesh
) const
{
    using namespace constant::mathematical;

    const clouds::spherical& sCloud = cloud<clouds::spherical>();
    const clouds::coupled& cCloud = cloud<clouds::coupled>();
    const clouds::sphericalCoupled& scCloud = cloud<clouds::sphericalCoupled>();

    const LagrangianSubVectorSubField U(cloud().U(model, subMesh));
    const LagrangianSubScalarSubField d(sCloud.d(model, subMesh));
    const LagrangianSubScalarField& v = sCloud.v(model, subMesh);
    const LagrangianSubScalarField& nuc = cCloud.nuc(model, subMesh);
    const LagrangianSubScalarField& Re = scCloud.Re(model, subMesh);
    const LagrangianSubVectorField& curlUc = cCloud.curlUc(model, subMesh);

    const LagrangianSubScalarField mcByMOrMc
    (
        isCloud<clouds::coupledToConstantDensityFluid>()
      ? v/cloud<clouds::coupledToConstantDensityFluid>().rhoByRhoc
      : v*cloud<clouds::coupledToFluid>().rhoc(model, subMesh)
    );

    const LagrangianSubScalarField Rew(mag(curlUc)*sqr(d)/nuc);
    const LagrangianSubScalarField beta(Rew/(Re + rootVSmall)/2);
    const LagrangianSubScalarField alpha(0.3314*sqrt(beta));
    const LagrangianSubScalarField Cld
    (
        neg(Re - 40)*6.46*((1 - alpha)*exp(-0.1*Re) + alpha)
      + pos(Re - 40)*6.46*0.0524*sqrt(beta*Re)
    );

    return
        LagrangianSubTensorField::New
        (
            subMesh.sub("L"),
            -mcByMOrMc*3/(twoPi*sqrt(Rew + rootVSmall))*Cld*(*curlUc)
        );
}


void Foam::Lagrangian::SaffmanMeiLift::addUSup
(
    const LagrangianSubVectorSubField& U,
    LagrangianEqn<vector>& eqn
) const
{
    const LagrangianSubTensorField& L = this->L(U.mesh());

    const LagrangianSubVectorField& Uc = cloud<clouds::carried>().Uc(U.mesh());

    if (eqn.isPsi(U))
    {
        eqn.Su += L & Uc;
        eqn -= Lagrangianm::Sp(L, U);
    }
    else
    {
        eqn += Lagrangianm::Sp(L, Uc);
        eqn.Su -= L & U;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Lagrangian::SaffmanMeiLift::SaffmanMeiLift
(
    const word& name,
    const LagrangianMesh& mesh,
    const dictionary& modelDict,
    const dictionary& stateDict
)
:
    LagrangianModel(name, mesh),
    cloudLagrangianModel(static_cast<const LagrangianModel&>(*this)),
    L
    (
        cloud().derivedField<tensor>
        (
            *this,
            &SaffmanMeiLift::calcL
        )
    )
{
    cloud<clouds::carried>().curlUc.psi();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::wordList Foam::Lagrangian::SaffmanMeiLift::addSupFields() const
{
    return wordList({cloud().U.name()});
}


bool Foam::Lagrangian::SaffmanMeiLift::addsSupToField
(
    const word& fieldName,
    const word& eqnFieldName
) const
{
    return
        fieldName == cloud().U.name()
     && (
            eqnFieldName == cloud().U.name()
         || eqnFieldName == cloud<clouds::coupled>().Uc.name()
        );
}


void Foam::Lagrangian::SaffmanMeiLift::addSup
(
    const LagrangianSubScalarField& deltaT,
    const LagrangianSubVectorSubField& U,
    LagrangianEqn<vector>& eqn
) const
{
    assertCloud<clouds::coupledToConstantDensityFluid>();

    addUSup(U, eqn);
}


void Foam::Lagrangian::SaffmanMeiLift::addSup
(
    const LagrangianSubScalarField& deltaT,
    const LagrangianSubScalarSubField& vOrM,
    const LagrangianSubVectorSubField& U,
    LagrangianEqn<vector>& eqn
) const
{
    assertCloud
    <
        clouds::coupledToConstantDensityFluid,
        clouds::coupledToFluid
    >();

    addUSup(U, eqn);
}


// ************************************************************************* //
