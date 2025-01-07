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

#include "cloudMassFlux.H"
#include "grouped.H"
#include "massive.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(cloudMassFlux, 0);
    addToRunTimeSelectionTable(functionObject, cloudMassFlux, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::LagrangianSubScalarSubField>
Foam::functionObjects::cloudMassFlux::q
(
    const LagrangianSubScalarSubField& fraction
) const
{
    return
        isCloud<clouds::grouped>()
      ? toSubField
        (
            cloud<clouds::grouped>().number(fraction.mesh())
           *cloud<clouds::massive>().m(fraction.mesh())
        )
      : toSubField
        (
            cloud<clouds::massive>().m(fraction.mesh())
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::cloudMassFlux::cloudMassFlux
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    cloudFlux(name, runTime, dict, clouds::massive::mName, dimMass)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::cloudMassFlux::~cloudMassFlux()
{}


// ************************************************************************* //
