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

#include "aerosolDrag.H"
#include "swarmCorrection.H"
#include "addToRunTimeSelectionTable.H"
#include "constants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(aerosolDrag, 0);
    addToRunTimeSelectionTable(dragModel, aerosolDrag, dictionary);
}
}

using Foam::constant::physicoChemical::k;
using Foam::constant::mathematical::pi;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::aerosolDrag::aerosolDrag
(
    const dictionary& dict,
    const phaseInterface& interface,
    const bool registerObject
)
:
    dispersedDragModel(dict, interface, registerObject),
    A1_(dict.lookupOrDefault<scalar>("A1", 2.514)),
    A2_(dict.lookupOrDefault<scalar>("A2", 0.8)),
    A3_(dict.lookupOrDefault<scalar>("A3", 0.55)),
    sigma_("sigma", dimLength, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModels::aerosolDrag::~aerosolDrag()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::dragModels::aerosolDrag::CdRe() const
{
    const volScalarField& T = interface_.continuous().thermo().T();
    const volScalarField& p = interface_.continuous().fluidThermo().p();
    tmp<volScalarField> td(interface_.dispersed().d());
    const volScalarField& d = td();

    const volScalarField lambda(k*T/(sqrt(2.0)*pi*p*sqr(sigma_)));

    return 24/(1 + lambda/d*(A1_ + A2_*exp(-A3_*d/lambda)));
}


// ************************************************************************* //
