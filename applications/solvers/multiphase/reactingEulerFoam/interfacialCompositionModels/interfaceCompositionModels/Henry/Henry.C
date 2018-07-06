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

#include "Henry.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::interfaceCompositionModels::Henry<Thermo, OtherThermo>::Henry
(
    const dictionary& dict,
    const phasePair& pair
)
:
    InterfaceCompositionModel<Thermo, OtherThermo>(dict, pair),
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
        dimensionedScalar("one", dimless, 1)
    )
{
    if (k_.size() != this->speciesNames_.size())
    {
        FatalErrorInFunction
            << "Differing number of species and solubilities"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::interfaceCompositionModels::Henry<Thermo, OtherThermo>::~Henry()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
void Foam::interfaceCompositionModels::Henry<Thermo, OtherThermo>::update
(
    const volScalarField& Tf
)
{
    YSolvent_ = scalar(1);

    forAllConstIter(hashedWordList, this->speciesNames_, iter)
    {
        YSolvent_ -= Yf(*iter, Tf);
    }
}


template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::interfaceCompositionModels::Henry<Thermo, OtherThermo>::Yf
(
    const word& speciesName,
    const volScalarField& Tf
) const
{
    if (this->speciesNames_.contains(speciesName))
    {
        const label index = this->speciesNames_[speciesName];

        return
            k_[index]
           *this->otherThermo_.composition().Y(speciesName)
           *this->otherThermo_.rhoThermo::rho()
           /this->thermo_.rhoThermo::rho();
    }
    else
    {
        return
            YSolvent_
           *this->thermo_.composition().Y(speciesName);
    }
}


template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::interfaceCompositionModels::Henry<Thermo, OtherThermo>::YfPrime
(
    const word& speciesName,
    const volScalarField& Tf
) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("YfPrime", this->pair_.name()),
                this->pair_.phase1().mesh().time().timeName(),
                this->pair_.phase1().mesh()
            ),
            this->pair_.phase1().mesh(),
            dimensionedScalar("zero", dimless/dimTemperature, 0)
        )
    );
}


// ************************************************************************* //
