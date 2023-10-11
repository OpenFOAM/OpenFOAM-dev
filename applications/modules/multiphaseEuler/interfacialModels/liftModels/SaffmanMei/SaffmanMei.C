/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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

#include "SaffmanMei.H"
#include "fvcCurl.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace liftModels
{
    defineTypeNameAndDebug(SaffmanMei, 0);
    addToRunTimeSelectionTable(liftModel, SaffmanMei, dictionary);
}
}

using Foam::constant::mathematical::twoPi;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::liftModels::SaffmanMei::SaffmanMei
(
    const dictionary& dict,
    const phaseInterface& interface
)
:
    dispersedLiftModel(dict, interface),
    residualRe_("residualRe", dimless, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::liftModels::SaffmanMei::~SaffmanMei()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::liftModels::SaffmanMei::Cl() const
{
    const volScalarField Re(max(interface_.Re(), residualRe_));
    const volScalarField Rew
    (
        mag(fvc::curl(interface_.continuous().U()))
       *sqr(interface_.dispersed().d())
       /(
            interface_.continuous().fluidThermo().nu()
          + dimensionedScalar(dimViscosity, small)
        )
    );

    const volScalarField Cld
    (
        neg0(Re - 40)*6.46
       *(
           (1 - 0.3314*sqrt(0.5*(Rew/Re)))*exp(-0.1*Re)
         + 0.3314*sqrt(0.5*(Rew/Re))
        )
      + pos(Re - 40)*6.46*0.0524*sqrt(0.5*(Rew/Re)*Re)
    );

    return 3/(twoPi*sqrt(Rew + small))*Cld;
}


// ************************************************************************* //
