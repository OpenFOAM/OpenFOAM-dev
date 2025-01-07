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

#include "coupledToIncompressibleFluid.H"
#include "viscosity.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace clouds
{
    defineTypeNameAndDebug(coupledToIncompressibleFluid, 0);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::clouds::coupledToIncompressibleFluid::getNucVf() const
{
    if (mesh_.mesh().foundObject<volScalarField>("nu"))
    {
        return mesh_.mesh().lookupObject<volScalarField>("nu");
    }

    if (mesh_.mesh().foundObject<viscosity>(physicalProperties::typeName))
    {
        return
            mesh_.mesh().lookupObject<viscosity>
            (
                physicalProperties::typeName
            ).nu();
    }

    FatalErrorInFunction
        << "Could not determine the carrier viscosity"
        << exit(FatalError);

    return tmp<volScalarField>(nullptr);
}


Foam::tmp<Foam::LagrangianSubScalarField>
Foam::clouds::coupledToIncompressibleFluid::calcNuc
(
    const LagrangianModelRef& model,
    const LagrangianSubMesh& subMesh
) const
{
    return nuc(model, subMesh);
}


// * * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * //

void Foam::clouds::coupledToIncompressibleFluid::updateNuc()
{
    if (tnucVf_.isTmp())
    {
        tnucVf_ = getNucVf();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::clouds::coupledToIncompressibleFluid::coupledToIncompressibleFluid
(
    const cloud& c
)
:
    coupled(c),
    mesh_(c.mesh()),
    physicalProperties_(c.mesh(), word::null),
    tnucVf_(getNucVf()),
    rhoByRhoc("rhoByRhoc", dimless, physicalProperties_),
    nuc(carrierField<scalar>(tnucVf_()))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::clouds::coupledToIncompressibleFluid::~coupledToIncompressibleFluid()
{}


// ************************************************************************* //
