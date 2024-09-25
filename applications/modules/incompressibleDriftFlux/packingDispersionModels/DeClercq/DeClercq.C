/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "DeClercq.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace packingDispersionModels
{
    defineTypeNameAndDebug(DeClercq, 0);
    addToRunTimeSelectionTable
    (
        packingDispersionModel,
        DeClercq,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::packingDispersionModels::DeClercq::DeClercq
(
    const dictionary& dict,
    const relativeVelocityModel& relativeVelocity
)
:
    packingDispersionModel(relativeVelocity),
    sigma0_("sigma0", sqr(dimVelocity), dict),
    beta_("beta", dimless, dict),
    alphaGel_("alphaGel", dimless, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::packingDispersionModels::DeClercq::~DeClercq()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::packingDispersionModels::DeClercq::sigmaPrime() const
{
    const volScalarField alphadmGel(mixture_.alphad() - alphaGel_);
    return pos(alphadmGel)*sigma0_/(beta_ + alphadmGel);
}


// ************************************************************************* //
