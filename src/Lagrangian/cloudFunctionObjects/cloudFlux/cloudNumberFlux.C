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

#include "cloudNumberFlux.H"
#include "grouped.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(cloudNumberFlux, 0);
    addToRunTimeSelectionTable(functionObject, cloudNumberFlux, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::LagrangianSubScalarSubField>
Foam::functionObjects::cloudNumberFlux::q
(
    const LagrangianSubScalarSubField& fraction
) const
{
    return
        isCloud<clouds::grouped>()
      ? cloud<clouds::grouped>().number(fraction.mesh())
      : toSubField<scalar, LagrangianSubMesh>
        (
            "1:" + Foam::name(fraction.mesh().group()),
            fraction.mesh(),
            dimensionedScalar(dimless, scalar(1))
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::cloudNumberFlux::cloudNumberFlux
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    cloudFlux(name, runTime, dict, "number", dimless)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::cloudNumberFlux::~cloudNumberFlux()
{}


// ************************************************************************* //
