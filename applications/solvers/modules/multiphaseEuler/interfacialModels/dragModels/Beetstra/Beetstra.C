/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2023 OpenFOAM Foundation
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

#include "Beetstra.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(Beetstra, 0);
    addToRunTimeSelectionTable(dragModel, Beetstra, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::Beetstra::Beetstra
(
    const dictionary& dict,
    const phaseInterface& interface,
    const bool registerObject
)
:
    dispersedDragModel(dict, interface, registerObject)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModels::Beetstra::~Beetstra()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::dragModels::Beetstra::CdRe() const
{
    const volScalarField alpha1
    (
        max(interface_.dispersed(), interface_.continuous().residualAlpha())
    );

    const volScalarField alpha2
    (
        max(1 - interface_.dispersed(), interface_.continuous().residualAlpha())
    );

    const volScalarField Res(max(alpha2*interface_.Re(), scalar(0)));

    const volScalarField F0
    (
        "F0",
        10*alpha1/sqr(alpha2) + sqr(alpha2)*(1 + 1.5*sqrt(alpha1))
    );

    const volScalarField F1
    (
        "F1",
        0.413*Res/(24*sqr(alpha2))*(1.0/alpha2
        + 3*alpha1*alpha2 + 8.4*pow(Res, -0.343))
        /(1 + pow(10.0, 3*alpha1)*pow(Res, -(1 + 4*alpha1)/2.0))
    );

    return 24*alpha2*(F0 + F1);
}


// ************************************************************************* //
