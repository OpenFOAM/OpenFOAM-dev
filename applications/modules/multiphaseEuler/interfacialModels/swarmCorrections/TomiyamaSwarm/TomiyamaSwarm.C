/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2023 OpenFOAM Foundation
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

#include "TomiyamaSwarm.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace swarmCorrections
{
    defineTypeNameAndDebug(TomiyamaSwarm, 0);
    addToRunTimeSelectionTable
    (
        swarmCorrection,
        TomiyamaSwarm,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::swarmCorrections::TomiyamaSwarm::TomiyamaSwarm
(
    const dictionary& dict,
    const phaseInterface& interface
)
:
    swarmCorrection(dict, interface),
    residualAlpha_
    (
        "residualAlpha",
        dimless,
        dict.lookupOrDefault<scalar>
        (
            "residualAlpha",
            interface_.dispersed().residualAlpha().value()
        )
    ),
    l_("l", dimless, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::swarmCorrections::TomiyamaSwarm::~TomiyamaSwarm()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::swarmCorrections::TomiyamaSwarm::Cs() const
{
    return pow(max(interface_.continuous(), residualAlpha_), 3 - 2*l_);
}


// ************************************************************************* //
