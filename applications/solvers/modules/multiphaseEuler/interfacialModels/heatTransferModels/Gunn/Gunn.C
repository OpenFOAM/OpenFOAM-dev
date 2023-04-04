/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2023 OpenFOAM Foundation
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
    const phaseInterface& interface,
    const bool registerObject
)
:
    heatTransferModel(dict, interface, registerObject),
    interface_
    (
        interface.modelCast<heatTransferModel, dispersedPhaseInterface>()
    )
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
        max(interface_.continuous(), interface_.continuous().residualAlpha())
    );

    const volScalarField sqrAlpha2(sqr(alpha2));

    const volScalarField Nu
    (
        (7 - 10*alpha2 + 5*sqrAlpha2)
       *(1 + 0.7*pow(interface_.Re(), 0.2)*cbrt(interface_.Pr()))
      + (1.33 - 2.4*alpha2 + 1.2*sqrAlpha2)
       *pow(interface_.Re(), 0.7)*cbrt(interface_.Pr())
    );

    return
        6*max(interface_.dispersed(), residualAlpha)
       *interface_.continuous().thermo().kappa()
       *Nu/sqr(interface_.dispersed().d());
}


// ************************************************************************* //
