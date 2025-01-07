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

#include "cloudBoundaryCollisionForce.H"
#include "grouped.H"
#include "massive.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(cloudBoundaryCollisionForce, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        cloudBoundaryCollisionForce,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::LagrangianSubScalarSubField>
Foam::functionObjects::cloudBoundaryCollisionForce::q
(
    const LagrangianSubScalarSubField& fraction,
    const label sign
) const
{
    const LagrangianSubMesh& subMesh = fraction.mesh();

    tmp<LagrangianSubScalarField> tp =
        sign
       *cloud<clouds::massive>().m(subMesh)
       *(
            cloud().U(subMesh)
          & fraction.mesh().nf<LagrangianSubVectorField>(fraction)
        );

    return
        toSubField
        (
            isCloud<clouds::grouped>()
          ? cloud<clouds::grouped>().number(subMesh)*tp
          : tp
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::cloudBoundaryCollisionForce::cloudBoundaryCollisionForce
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    cloudBoundaryCollisionFlux(name, runTime, dict, "F", dimForce)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::cloudBoundaryCollisionForce::
~cloudBoundaryCollisionForce()
{}


// ************************************************************************* //
