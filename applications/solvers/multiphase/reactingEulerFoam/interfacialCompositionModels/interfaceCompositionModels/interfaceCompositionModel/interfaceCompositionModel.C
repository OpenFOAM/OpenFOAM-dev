/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2019 OpenFOAM Foundation
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
#include "phasePair.H"
#include "rhoReactionThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(interfaceCompositionModel, 0);
    defineRunTimeSelectionTable(interfaceCompositionModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceCompositionModel::interfaceCompositionModel
(
    const dictionary& dict,
    const phasePair& pair
)
:
    pair_(pair),
    species_(dict.lookup("species")),
    Le_("Le", dimless, dict),
    thermo_
    (
        pair.phase1().mesh().lookupObject<rhoReactionThermo>
        (
            IOobject::groupName(basicThermo::dictName, pair.phase1().name())
        )
    ),
    otherThermo_
    (
        pair.phase2().mesh().lookupObject<rhoThermo>
        (
            IOobject::groupName(basicThermo::dictName, pair.phase2().name())
        )
    )
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
        IOobject::groupName("D", pair_.name()),
        composition().alphah(speciei, p, T)
       /composition().rho(speciei, p, T)
       /Le_
    );
}


Foam::tmp<Foam::volScalarField> Foam::interfaceCompositionModel::L
(
    const word& speciesName,
    const volScalarField& Tf
) const
{
    const label speciei = composition().species()[speciesName];
    const volScalarField& p(thermo_.p());
    volScalarField Ha(composition().Ha(speciei, p, Tf));

    const volScalarField& otherP(otherThermo_.p());
    tmp<volScalarField> otherHa(nullptr);
    if (otherHasComposition())
    {
        const label otherSpeciei = otherComposition().species()[speciesName];
        otherHa = otherComposition().Ha(otherSpeciei, otherP, Tf);
    }
    else
    {
        otherHa = otherThermo_.ha(otherP, Tf);
    }

    return
        volScalarField::New
        (
            IOobject::groupName("L", pair_.name()),
            otherHa - Ha
        );
}


void Foam::interfaceCompositionModel::addDmdtL
(
    const volScalarField& K,
    const volScalarField& Tf,
    volScalarField& dmdtL,
    volScalarField& dmdtLPrime
) const
{
    forAllConstIter(hashedWordList, species_, iter)
    {
        const volScalarField rhoKDL(thermo_.rho()*K*D(*iter)*L(*iter, Tf));

        dmdtL += rhoKDL*dY(*iter, Tf);
        dmdtLPrime += rhoKDL*YfPrime(*iter, Tf);
    }
}


// ************************************************************************* //
