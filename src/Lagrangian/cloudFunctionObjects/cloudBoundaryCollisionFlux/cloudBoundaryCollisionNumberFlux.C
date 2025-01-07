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

#include "cloudBoundaryCollisionNumberFlux.H"
#include "grouped.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(cloudBoundaryCollisionNumberFlux, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        cloudBoundaryCollisionNumberFlux,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::LagrangianSubScalarSubField>
Foam::functionObjects::cloudBoundaryCollisionNumberFlux::q
(
    const LagrangianSubScalarSubField& fraction,
    const label sign
) const
{
    const LagrangianSubMesh& subMesh = fraction.mesh();

    if (sign > 0 && isCloud<clouds::grouped>())
    {
        return cloud<clouds::grouped>().number(subMesh);
    }

    return
        toSubField<scalar, LagrangianSubMesh>
        (
            Foam::name(scalar(sign > 0)) + ":" + Foam::name(subMesh.group()),
            subMesh,
            dimensionedScalar(dimless, scalar(sign > 0))
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::cloudBoundaryCollisionNumberFlux::
cloudBoundaryCollisionNumberFlux
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    cloudBoundaryCollisionFlux(name, runTime, dict, "phiNumber", inv(dimTime))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::cloudBoundaryCollisionNumberFlux::
~cloudBoundaryCollisionNumberFlux()
{}


// ************************************************************************* //
