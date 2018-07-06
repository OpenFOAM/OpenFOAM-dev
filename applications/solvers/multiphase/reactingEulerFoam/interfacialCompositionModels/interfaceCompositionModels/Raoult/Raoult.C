/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::interfaceCompositionModels::Raoult<Thermo, OtherThermo>::Raoult
(
    const dictionary& dict,
    const phasePair& pair
)
:
    InterfaceCompositionModel<Thermo, OtherThermo>(dict, pair),
    YNonVapour_
    (
        IOobject
        (
            IOobject::groupName("YNonVapour", pair.name()),
            pair.phase1().mesh().time().timeName(),
            pair.phase1().mesh()
        ),
        pair.phase1().mesh(),
        dimensionedScalar("one", dimless, 1)
    ),
    YNonVapourPrime_
    (
        IOobject
        (
            IOobject::groupName("YNonVapourPrime", pair.name()),
            pair.phase1().mesh().time().timeName(),
            pair.phase1().mesh()
        ),
        pair.phase1().mesh(),
        dimensionedScalar("zero", dimless/dimTemperature, 0)
    )
{
    forAllConstIter(hashedWordList, this->speciesNames_, iter)
    {
        speciesModels_.insert
        (
            *iter,
            autoPtr<interfaceCompositionModel>
            (
                interfaceCompositionModel::New
                (
                    dict.subDict(*iter),
                    pair
                )
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::interfaceCompositionModels::Raoult<Thermo, OtherThermo>::~Raoult()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
void Foam::interfaceCompositionModels::Raoult<Thermo, OtherThermo>::update
(
    const volScalarField& Tf
)
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
            this->otherThermo_.composition().Y(iter.key())
           *iter()->Yf(iter.key(), Tf);

        YNonVapourPrime_ -=
            this->otherThermo_.composition().Y(iter.key())
           *iter()->YfPrime(iter.key(), Tf);
    }
}


template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::interfaceCompositionModels::Raoult<Thermo, OtherThermo>::Yf
(
    const word& speciesName,
    const volScalarField& Tf
) const
{
    if (this->speciesNames_.contains(speciesName))
    {
        return
             this->otherThermo_.composition().Y(speciesName)
            *speciesModels_[speciesName]->Yf(speciesName, Tf);
    }
    else
    {
        return
             this->thermo_.composition().Y(speciesName)
            *YNonVapour_;
    }
}


template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::interfaceCompositionModels::Raoult<Thermo, OtherThermo>::YfPrime
(
    const word& speciesName,
    const volScalarField& Tf
) const
{
    if (this->speciesNames_.contains(speciesName))
    {
        return
             this->otherThermo_.composition().Y(speciesName)
            *speciesModels_[speciesName]->YfPrime(speciesName, Tf);
    }
    else
    {
        return
            this->otherThermo_.composition().Y(speciesName)
           *YNonVapourPrime_;
    }
}


// ************************************************************************* //
