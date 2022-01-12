/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2022 OpenFOAM Foundation
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

#include "interfaceCompositionModel.H"
#include "phaseModel.H"
#include "phaseSystem.H"
#include "rhoReactionThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(interfaceCompositionModel, 0);
    defineSidedInterfacialModelTypeNameAndDebug(interfaceCompositionModel, 0);
    defineRunTimeSelectionTable(interfaceCompositionModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceCompositionModel::interfaceCompositionModel
(
    const dictionary& dict,
    const phaseInterface& interface
)
:
    interface_
    (
        interface.modelCast<interfaceCompositionModel, sidedPhaseInterface>()
    ),
    species_(dict.lookup("species")),
    Le_("Le", dimless, dict),
    thermo_(refCast<const rhoReactionThermo>(interface_.phase().thermo())),
    otherThermo_(interface_.otherPhase().thermo())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interfaceCompositionModel::~interfaceCompositionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::interfaceCompositionModel::dY
(
    const word& speciesName,
    const volScalarField& Tf
) const
{
    const label speciei = composition().species()[speciesName];

    return Yf(speciesName, Tf) - composition().Y()[speciei];
}


Foam::tmp<Foam::volScalarField> Foam::interfaceCompositionModel::dYfPrime
(
    const word& speciesName,
    const volScalarField& Tf
) const
{
    return YfPrime(speciesName, Tf);
}


Foam::tmp<Foam::volScalarField> Foam::interfaceCompositionModel::D
(
    const word& speciesName
) const
{
    const label speciei = composition().species()[speciesName];
    const volScalarField& p(thermo_.p());
    const volScalarField& T(thermo_.T());

    return volScalarField::New
    (
        IOobject::groupName("D" + speciesName, interface_.name()),
        composition().kappa(speciei, p, T)
       /composition().Cp(speciei, p, T)
       /composition().rho(speciei, p, T)
       /Le_
    );
}


// ************************************************************************* //
