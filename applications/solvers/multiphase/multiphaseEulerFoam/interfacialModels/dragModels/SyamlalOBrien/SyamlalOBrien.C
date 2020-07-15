/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "SyamlalOBrien.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(SyamlalOBrien, 0);
    addToRunTimeSelectionTable(dragModel, SyamlalOBrien, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::SyamlalOBrien::SyamlalOBrien
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject
)
:
    dragModel(dict, pair, registerObject)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModels::SyamlalOBrien::~SyamlalOBrien()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::dragModels::SyamlalOBrien::CdRe() const
{
    const volScalarField alpha2
    (
        max(1 - pair_.dispersed(), pair_.continuous().residualAlpha())
    );

    const volScalarField A(pow(alpha2, 4.14));
    const volScalarField B
    (
        neg(alpha2 - 0.85)*(0.8*pow(alpha2, 1.28))
      + pos0(alpha2 - 0.85)*(pow(alpha2, 2.65))
    );
    const volScalarField Re(pair_.Re());
    const volScalarField Vr
    (
        0.5
       *(
            A - 0.06*Re + sqrt(sqr(0.06*Re) + 0.12*Re*(2*B - A) + sqr(A))
        )
    );
    volScalarField CdsRe(sqr(0.63*sqrt(Re) + 4.8*sqrt(Vr)));

    return
        CdsRe
       *max(pair_.continuous(), pair_.continuous().residualAlpha())
       /sqr(Vr);
}


// ************************************************************************* //
