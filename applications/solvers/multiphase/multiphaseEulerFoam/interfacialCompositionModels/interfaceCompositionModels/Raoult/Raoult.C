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

#include "Raoult.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace interfaceCompositionModels
{
    defineTypeNameAndDebug(Raoult, 0);
    addToRunTimeSelectionTable(interfaceCompositionModel, Raoult, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceCompositionModels::Raoult::Raoult
(
    const dictionary& dict,
    const phaseInterface& interface
)
:
    interfaceCompositionModel(dict, interface),
    YNonVapour_
    (
        IOobject
        (
            IOobject::groupName("YNonVapour", this->interface().name()),
            interface.mesh().time().timeName(),
            interface.mesh()
        ),
        interface.mesh(),
        dimensionedScalar(dimless, 1)
    ),
    YNonVapourPrime_
    (
        IOobject
        (
            IOobject::groupName("YNonVapourPrime", this->interface().name()),
            interface.mesh().time().timeName(),
            interface.mesh()
        ),
        interface.mesh(),
        dimensionedScalar(dimless/dimTemperature, 0)
    )
{
    forAllConstIter(hashedWordList, species(), iter)
    {
        speciesModels_.insert
        (
            *iter,
            autoPtr<interfaceCompositionModel>
            (
                interfaceCompositionModel::New
                (
                    dict.subDict(*iter),
                    interface
                )
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interfaceCompositionModels::Raoult::~Raoult()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::interfaceCompositionModels::Raoult::update(const volScalarField& Tf)
{
    YNonVapour_ = scalar(1);

    forAllIter
    (
        HashTable<autoPtr<interfaceCompositionModel>>,
        speciesModels_,
        iter
    )
    {
        iter()->update(Tf);

        YNonVapour_ -=
            otherComposition().Y(iter.key())
           *iter()->Yf(iter.key(), Tf);

        YNonVapourPrime_ -=
            otherComposition().Y(iter.key())
           *iter()->YfPrime(iter.key(), Tf);
    }
}


Foam::tmp<Foam::volScalarField> Foam::interfaceCompositionModels::Raoult::Yf
(
    const word& speciesName,
    const volScalarField& Tf
) const
{
    if (species().found(speciesName))
    {
        return
            otherComposition().Y(speciesName)
           *speciesModels_[speciesName]->Yf(speciesName, Tf);
    }
    else
    {
        return composition().Y(speciesName)*YNonVapour_;
    }
}


Foam::tmp<Foam::volScalarField>
Foam::interfaceCompositionModels::Raoult::YfPrime
(
    const word& speciesName,
    const volScalarField& Tf
) const
{
    if (species().found(speciesName))
    {
        return
            otherComposition().Y(speciesName)
           *speciesModels_[speciesName]->YfPrime(speciesName, Tf);
    }
    else
    {
        return composition().Y(speciesName)*YNonVapourPrime_;
    }
}


// ************************************************************************* //
