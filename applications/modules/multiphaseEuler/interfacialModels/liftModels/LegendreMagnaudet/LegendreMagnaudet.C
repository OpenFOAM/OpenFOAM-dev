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

#include "LegendreMagnaudet.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace liftModels
{
    defineTypeNameAndDebug(LegendreMagnaudet, 0);
    addToRunTimeSelectionTable(liftModel, LegendreMagnaudet, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::liftModels::LegendreMagnaudet::LegendreMagnaudet
(
    const dictionary& dict,
    const phaseInterface& interface
)
:
    dispersedLiftModel(dict, interface),
    residualRe_("residualRe", dimless, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::liftModels::LegendreMagnaudet::~LegendreMagnaudet()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::liftModels::LegendreMagnaudet::Cl() const
{
    const volScalarField Re(max(interface_.Re(), residualRe_));

    const volScalarField Sr
    (
        sqr(interface_.dispersed().d())
       /(Re*interface_.continuous().fluidThermo().nu())
       *mag(fvc::grad(interface_.continuous().U()))
    );

    const volScalarField ClLowSqr
    (
        sqr(6*2.255)*sqr(Sr)
       /(pow4(constant::mathematical::pi)*Re*pow3(Sr + 0.2*Re))
    );

    const volScalarField ClHighSqr(sqr(0.5*(Re + 16)/(Re + 29)));

    return sqrt(ClLowSqr + ClHighSqr);
}


// ************************************************************************* //
