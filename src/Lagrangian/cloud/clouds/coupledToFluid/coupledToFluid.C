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

#include "coupledToFluid.H"
#include "fluidThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace clouds
{
    defineTypeNameAndDebug(coupledToFluid, 0);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::clouds::coupledToFluid::getRhocVf() const
{
    if (mesh_.mesh().foundObject<volScalarField>("rho"))
    {
        return mesh_.mesh().lookupObject<volScalarField>("rho");
    }

    if (mesh_.mesh().foundObject<basicThermo>(physicalProperties::typeName))
    {
        return
            mesh_.mesh().lookupObject<basicThermo>
            (
                physicalProperties::typeName
            ).rho();
    }

    FatalErrorInFunction
        << "Could not determine the carrier density"
        << exit(FatalError);

    return tmp<volScalarField>(nullptr);
}


Foam::tmp<Foam::LagrangianSubScalarField>
Foam::clouds::coupledToFluid::calcNuc
(
    const LagrangianModelRef& model,
    const LagrangianSubMesh& subMesh
) const
{
    return muc(model, subMesh)/rhoc(model, subMesh);
}


// * * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * //

void Foam::clouds::coupledToFluid::updateRhoc()
{
    if (trhocVf_.isTmp())
    {
        trhocVf_ = getRhocVf();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::clouds::coupledToFluid::coupledToFluid(const cloud& c)
:
    coupled(c),
    mesh_(c.mesh()),
    trhocVf_(getRhocVf()),
    rhoc(carrierField<scalar>(trhocVf_())),
    muc
    (
        carrierField<scalar>
        (
            mesh_
           .mesh()
           .lookupObject<fluidThermo>(physicalProperties::typeName)
           .mu()
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::clouds::coupledToFluid::~coupledToFluid()
{}


// ************************************************************************* //
