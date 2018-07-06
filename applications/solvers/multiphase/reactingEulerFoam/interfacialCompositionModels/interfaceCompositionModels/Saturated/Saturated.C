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

#include "Saturated.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::interfaceCompositionModels::Saturated<Thermo, OtherThermo>::
wRatioByP() const
{
    const dimensionedScalar Wi
    (
        "W",
        dimMass/dimMoles,
        this->thermo_.composition().W(saturatedIndex_)
    );

    return Wi/this->thermo_.W()/this->thermo_.p();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::interfaceCompositionModels::Saturated<Thermo, OtherThermo>::Saturated
(
    const dictionary& dict,
    const phasePair& pair
)
:
    InterfaceCompositionModel<Thermo, OtherThermo>(dict, pair),
    saturatedName_(this->speciesNames_[0]),
    saturatedIndex_
    (
        this->thermo_.composition().species()[saturatedName_]
    ),
    saturationModel_
    (
        saturationModel::New
        (
            dict.subDict("saturationPressure"),
            pair.phase1().mesh()
        )
    )
{
    if (this->speciesNames_.size() != 1)
    {
        FatalErrorInFunction
            << "Saturated model is suitable for one species only."
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::interfaceCompositionModels::Saturated<Thermo, OtherThermo>::~Saturated()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
void
Foam::interfaceCompositionModels::Saturated<Thermo, OtherThermo>::update
(
    const volScalarField& Tf
)
{}


template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::interfaceCompositionModels::Saturated<Thermo, OtherThermo>::Yf
(
    const word& speciesName,
    const volScalarField& Tf
) const
{
    if (saturatedName_ == speciesName)
    {
        return wRatioByP()*saturationModel_->pSat(Tf);
    }
    else
    {
        const label speciesIndex
        (
            this->thermo_.composition().species()[speciesName]
        );

        return
            this->thermo_.Y()[speciesIndex]
           *(scalar(1) - wRatioByP()*saturationModel_->pSat(Tf))
           /max(scalar(1) - this->thermo_.Y()[saturatedIndex_], small);
    }
}


template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::interfaceCompositionModels::Saturated<Thermo, OtherThermo>::YfPrime
(
    const word& speciesName,
    const volScalarField& Tf
) const
{
    if (saturatedName_ == speciesName)
    {
        return wRatioByP()*saturationModel_->pSatPrime(Tf);
    }
    else
    {
        const label speciesIndex
        (
            this->thermo_.composition().species()[speciesName]
        );

        return
          - this->thermo_.Y()[speciesIndex]
           *wRatioByP()*saturationModel_->pSatPrime(Tf)
           /max(scalar(1) - this->thermo_.Y()[saturatedIndex_], small);
    }
}


// ************************************************************************* //
