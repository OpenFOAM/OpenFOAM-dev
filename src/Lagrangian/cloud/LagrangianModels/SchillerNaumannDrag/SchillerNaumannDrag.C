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

#include "SchillerNaumannDrag.H"
#include "addToRunTimeSelectionTable.H"
#include "coupledToConstantDensityFluid.H"
#include "coupledToFluid.H"
#include "sphericalCoupled.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace Lagrangian
{
    defineTypeNameAndDebug(SchillerNaumannDrag, 0);
    addToRunTimeSelectionTable
    (
        LagrangianModel,
        SchillerNaumannDrag,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::LagrangianSubScalarField>
Foam::Lagrangian::SchillerNaumannDrag::calcD
(
    const LagrangianModelRef& model,
    const LagrangianSubMesh& subMesh
) const
{
    const clouds::spherical& sCloud = cloud<clouds::spherical>();
    const clouds::sphericalCoupled& scCloud = cloud<clouds::sphericalCoupled>();

    const LagrangianSubScalarField& Re = scCloud.Re(model, subMesh);
    const LagrangianSubScalarSubField d(sCloud.d(model, subMesh));

    assertCloud
    <
        clouds::coupledToConstantDensityFluid,
        clouds::coupledToFluid
    >();

    tmp<LagrangianSubScalarField> tmucByRhoOrMuc =
        isCloud<clouds::coupledToConstantDensityFluid>()
      ? (
            cloud<clouds::coupledToConstantDensityFluid>().nuc(model, subMesh)
           /cloud<clouds::coupledToConstantDensityFluid>().rhoByRhoc
        )
      : tmp<LagrangianSubScalarField>
        (
            cloud<clouds::coupledToFluid>().muc(model, subMesh)
        );

    return
        LagrangianSubScalarField::New
        (
            subMesh.sub("D"),
            CdRe(Re)*(constant::mathematical::pi/8)*d*tmucByRhoOrMuc
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Lagrangian::SchillerNaumannDrag::SchillerNaumannDrag
(
    const word& name,
    const LagrangianMesh& mesh,
    const dictionary& modelDict,
    const dictionary& stateDict
)
:
    drag(name, mesh, modelDict, stateDict)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::LagrangianSubScalarField>
Foam::Lagrangian::SchillerNaumannDrag::CdRe
(
    const LagrangianSubScalarField& Re
)
{
    return
        neg(Re - 1000)*24*(1 + 0.15*pow(Re, 0.687))
      + pos0(Re - 1000)*0.44*Re;
}


// ************************************************************************* //
