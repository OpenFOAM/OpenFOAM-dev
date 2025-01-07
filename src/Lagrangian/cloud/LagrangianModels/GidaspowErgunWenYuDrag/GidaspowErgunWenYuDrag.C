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

#include "GidaspowErgunWenYuDrag.H"
#include "SchillerNaumannDrag.H"
#include "addToRunTimeSelectionTable.H"
#include "coupledToIncompressibleFluid.H"
#include "coupledToFluid.H"
#include "sphericalCoupled.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace Lagrangian
{
    defineTypeNameAndDebug(GidaspowErgunWenYuDrag, 0);
    addToRunTimeSelectionTable
    (
        LagrangianModel,
        GidaspowErgunWenYuDrag,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::LagrangianSubScalarField>
Foam::Lagrangian::GidaspowErgunWenYuDrag::calcD
(
    const LagrangianModelRef& model,
    const LagrangianSubMesh& subMesh
) const
{
    const clouds::spherical& sCloud = cloud<clouds::spherical>();
    const clouds::sphericalCoupled& scCloud = cloud<clouds::sphericalCoupled>();

    const LagrangianSubScalarField Re = scCloud.Re(model, subMesh);
    const LagrangianSubScalarSubField d(sCloud.d(model, subMesh));

    const LagrangianSubScalarField alpha(min(sCloud.alpha(subMesh), alphaMax_));
    const LagrangianSubScalarField alphac(1 - alpha);

    const LagrangianSubScalarField CdRe
    (
        // Use Wen-Yu at low particulate fractions (< 20%) ...
        pos0(alphac - 0.8)
       *SchillerNaumannDrag::CdRe(alphac*Re)
       *pow(alphac, -2.65)

        // ... and Ergun at high particulate fractions (> 20%)
      + neg(alphac - 0.8)*(4.0/3.0)*(150*alpha/alphac + 1.75*Re)
    );

    assertCloud
    <
        clouds::coupledToIncompressibleFluid,
        clouds::coupledToFluid
    >();

    tmp<LagrangianSubScalarField> tmucByRhoOrMuc =
        isCloud<clouds::coupledToIncompressibleFluid>()
      ? (
            cloud<clouds::coupledToIncompressibleFluid>().nuc(model, subMesh)
           /cloud<clouds::coupledToIncompressibleFluid>().rhoByRhoc
        )
      : tmp<LagrangianSubScalarField>
        (
            cloud<clouds::coupledToFluid>().muc(model, subMesh)
        );

    return
        LagrangianSubScalarField::New
        (
            "D:" + Foam::name(subMesh.group()),
            CdRe*(constant::mathematical::pi/8)*d*tmucByRhoOrMuc
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Lagrangian::GidaspowErgunWenYuDrag::GidaspowErgunWenYuDrag
(
    const word& name,
    const LagrangianMesh& mesh,
    const dictionary& modelDict,
    const dictionary& stateDict
)
:
    drag(name, mesh, modelDict, stateDict),
    alphaMax_(modelDict.lookup<scalar>("alphaMax", unitFraction))
{}


// ************************************************************************* //
