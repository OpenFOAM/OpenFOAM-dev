/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2024 OpenFOAM Foundation
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

#include "MichaelsBolger.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace relativeVelocityModels
{
    defineTypeNameAndDebug(MichaelsBolger, 0);
    addToRunTimeSelectionTable
    (
        relativeVelocityModel,
        MichaelsBolger,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativeVelocityModels::MichaelsBolger::MichaelsBolger
(
    const dictionary& dict,
    const incompressibleDriftFluxMixture& mixture,
    const uniformDimensionedVectorField& g
)
:
    relativeVelocityModel(dict, mixture, g),
    a0_("a0", dimless, dict),
    a1_("a1", dimless, dict),
    Vc_("Vc", dimTime, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::relativeVelocityModels::MichaelsBolger::~MichaelsBolger()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::relativeVelocityModels::MichaelsBolger::UdmCoeff() const
{
    return
        (mixture_.rhoc()/mixture_.rho())*Vc_
       *(
           a0_
         + pow(max(1 - mixture_.alphad()/mixture().alphaMax(), scalar(0)), a1_)
        );
}


// ************************************************************************* //
