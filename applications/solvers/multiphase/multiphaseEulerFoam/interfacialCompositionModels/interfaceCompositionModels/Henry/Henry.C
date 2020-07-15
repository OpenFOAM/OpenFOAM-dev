/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2020 OpenFOAM Foundation
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

#include "Henry.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace interfaceCompositionModels
{
    defineTypeNameAndDebug(Henry, 0);
    addToRunTimeSelectionTable(interfaceCompositionModel, Henry, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceCompositionModels::Henry::Henry
(
    const dictionary& dict,
    const phasePair& pair
)
:
    interfaceCompositionModel(dict, pair),
    k_(dict.lookup("k")),
    YSolvent_
    (
        IOobject
        (
            IOobject::groupName("YSolvent", pair.name()),
            pair.phase1().mesh().time().timeName(),
            pair.phase1().mesh()
        ),
        pair.phase1().mesh(),
        dimensionedScalar(dimless, 1)
    )
{
    if (k_.size() != species().size())
    {
        FatalErrorInFunction
            << "Differing number of species and solubilities"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interfaceCompositionModels::Henry::~Henry()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::interfaceCompositionModels::Henry::update(const volScalarField& Tf)
{
    YSolvent_ = scalar(1);

    forAllConstIter(hashedWordList, species(), iter)
    {
        YSolvent_ -= Yf(*iter, Tf);
    }
}


Foam::tmp<Foam::volScalarField> Foam::interfaceCompositionModels::Henry::Yf
(
    const word& speciesName,
    const volScalarField& Tf
) const
{
    if (species().found(speciesName))
    {
        const label index = species()[speciesName];

        return
            k_[index]
           *otherComposition().Y(speciesName)
           *otherThermo().rho()
           /thermo().rho();
    }
    else
    {
        return YSolvent_*composition().Y(speciesName);
    }
}


Foam::tmp<Foam::volScalarField> Foam::interfaceCompositionModels::Henry::YfPrime
(
    const word& speciesName,
    const volScalarField& Tf
) const
{
    return volScalarField::New
    (
        IOobject::groupName("YfPrime", pair().name()),
        pair().phase1().mesh(),
        dimensionedScalar(dimless/dimTemperature, 0)
    );
}


// ************************************************************************* //
