/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "RanzMarshallHeatTransfer.H"
#include "addToRunTimeSelectionTable.H"
#include "coupledToThermalFluid.H"
#include "sphericalCoupled.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace Lagrangian
{
    defineTypeNameAndDebug(RanzMarshallHeatTransfer, 0);
    addToRunTimeSelectionTable
    (
        LagrangianModel,
        RanzMarshallHeatTransfer,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::LagrangianSubScalarField>
Foam::Lagrangian::RanzMarshallHeatTransfer::calcH
(
    const LagrangianModelRef& model,
    const LagrangianSubMesh& subMesh
) const
{
    const clouds::spherical& sCloud = cloud<clouds::spherical>();
    const clouds::sphericalCoupled& scCloud = cloud<clouds::sphericalCoupled>();
    const clouds::coupledToThermalFluid& cctfCloud =
        cloud<clouds::coupledToThermalFluid>();

    tmp<LagrangianSubScalarSubField> td = sCloud.d(model, subMesh);
    const LagrangianSubScalarSubField& d = td();
    const LagrangianSubScalarField& a = sCloud.a(model, subMesh);
    const LagrangianSubScalarField& Re = scCloud.Re(model, subMesh);
    const LagrangianSubScalarField& kappac = cctfCloud.kappac(model, subMesh);
    const LagrangianSubScalarField& Prc = cctfCloud.Prc(model, subMesh);

    tmp<LagrangianSubScalarField> tNu(2 + 0.6*sqrt(Re)*cbrt(Prc));

    return a*tNu*kappac/d;
}


// ************************************************************************* //
