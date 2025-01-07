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

#include "cloudVolumeFlux.H"
#include "grouped.H"
#include "shaped.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(cloudVolumeFlux, 0);
    addToRunTimeSelectionTable(functionObject, cloudVolumeFlux, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::LagrangianSubScalarSubField>
Foam::functionObjects::cloudVolumeFlux::q
(
    const LagrangianSubScalarSubField& fraction
) const
{
    tmp<LagrangianSubScalarField> tv =
        cloud<clouds::shaped>().v(fraction.mesh());

    return
        toSubField
        (
            isCloud<clouds::grouped>()
          ? cloud<clouds::grouped>().number(fraction.mesh())*tv
          : tv
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::cloudVolumeFlux::cloudVolumeFlux
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    cloudFlux(name, runTime, dict, clouds::shaped::vName, dimVolume)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::cloudVolumeFlux::~cloudVolumeFlux()
{}


// ************************************************************************* //
