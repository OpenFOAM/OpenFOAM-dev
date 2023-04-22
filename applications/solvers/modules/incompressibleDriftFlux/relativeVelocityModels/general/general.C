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

#include "general.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace relativeVelocityModels
{
    defineTypeNameAndDebug(general, 0);
    addToRunTimeSelectionTable(relativeVelocityModel, general, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::relativeVelocityModels::general::general
(
    const dictionary& dict,
    const incompressibleDriftFluxMixture& mixture,
    const uniformDimensionedVectorField& g
)
:
    relativeVelocityModel(dict, mixture, g),
    a_("a", dimless, dict),
    a1_("a1", dimless, dict),
    Vc_("Vc", dimTime, dict),
    residualAlpha_("residualAlpha", dimless, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::relativeVelocityModels::general::~general()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::relativeVelocityModels::general::correct()
{
    const volScalarField& alphad = mixture_.alphad();

    Udm_ =
        (mixture_.rhoc()/mixture_.rho())*Vc_*acceleration()
       *(
            exp(-a_*max(alphad - residualAlpha_, scalar(0)))
          - exp(-a1_*max(alphad - residualAlpha_, scalar(0)))
        );
}


// ************************************************************************* //
