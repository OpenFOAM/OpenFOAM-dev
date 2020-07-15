/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2020 OpenFOAM Foundation
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

#include "Gunn.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferModels
{
    defineTypeNameAndDebug(Gunn, 0);
    addToRunTimeSelectionTable(heatTransferModel, Gunn, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatTransferModels::Gunn::Gunn
(
    const dictionary& dict,
    const phasePair& pair
)
:
    heatTransferModel(dict, pair)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::heatTransferModels::Gunn::~Gunn()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::heatTransferModels::Gunn::K(const scalar residualAlpha) const
{
    const volScalarField alpha2
    (
        max(1 - pair_.dispersed(), pair_.continuous().residualAlpha())
    );

    const volScalarField sqrAlpha2(sqr(alpha2));

    const volScalarField Nu
    (
        (7 - 10*alpha2 + 5*sqrAlpha2)
       *(1 + 0.7*pow(pair_.Re(), 0.2)*cbrt(pair_.Pr()))
      + (1.33 - 2.4*alpha2 + 1.2*sqrAlpha2)
       *pow(pair_.Re(), 0.7)*cbrt(pair_.Pr())
    );

    return
        6*max(pair_.dispersed(), residualAlpha)
       *pair_.continuous().thermo().kappa()
       *Nu/sqr(pair_.dispersed().d());
}


// ************************************************************************* //
