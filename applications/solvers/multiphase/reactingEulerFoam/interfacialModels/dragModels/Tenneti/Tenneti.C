/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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

#include "Tenneti.H"
#include "phasePair.H"
#include "SchillerNaumann.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(Tenneti, 0);
    addToRunTimeSelectionTable(dragModel, Tenneti, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::Tenneti::Tenneti
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject
)
:
    dragModel(dict, pair, registerObject),
    residualRe_("residualRe", dimless, dict.lookup("residualRe"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModels::Tenneti::~Tenneti()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::dragModels::Tenneti::CdRe() const
{
    volScalarField alpha1
    (
        max(pair_.dispersed(), pair_.continuous().residualAlpha())
    );

    volScalarField alpha2
    (
        max(scalar(1) - pair_.dispersed(), pair_.continuous().residualAlpha())
    );

    volScalarField Res(alpha2*pair_.Re());

    volScalarField CdReIsolated
    (
        neg(Res - 1000)*24*(1 + 0.15*pow(Res, 0.687))
      + pos0(Res - 1000)*0.44*max(Res, residualRe_)
    );

    volScalarField F0
    (
        5.81*alpha1/pow3(alpha2) + 0.48*pow(alpha1, 1.0/3.0)/pow4(alpha2)
    );

    volScalarField F1
    (
        pow3(alpha1)*Res*(0.95 + 0.61*pow3(alpha1)/sqr(alpha2))
    );

    // Tenneti et al. correlation includes the mean pressure drag.
    // This was removed here by multiplying F by alpha2 for consistency with
    // the formulation used in OpenFOAM
    return
        CdReIsolated + 24*sqr(alpha2)*(F0 + F1);
}


// ************************************************************************* //
