/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2023 OpenFOAM Foundation
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

#include "saturated.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace interfaceCompositionModels
{
    defineTypeNameAndDebug(saturated, 0);
    addToRunTimeSelectionTable
    (
        interfaceCompositionModel,
        saturated,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::interfaceCompositionModels::saturated::wRatioByP() const
{
    const dimensionedScalar Wi
    (
        "W",
        dimMass/dimMoles,
        thermo().Wi(saturatedIndex_)
    );

    return Wi/thermo().W()/thermo().p();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceCompositionModels::saturated::saturated
(
    const dictionary& dict,
    const phaseInterface& interface
)
:
    interfaceCompositionModel(dict, interface),
    saturatedName_(species()[0]),
    saturatedIndex_(thermo().species()[saturatedName_]),
    saturationModel_(saturationPressureModel::New("pSat", dict))
{
    if (species().size() != 1)
    {
        FatalErrorInFunction
            << "saturated model is suitable for one species only."
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interfaceCompositionModels::saturated::~saturated()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::interfaceCompositionModels::saturated::update
(
    const volScalarField& Tf
)
{}


Foam::tmp<Foam::volScalarField> Foam::interfaceCompositionModels::saturated::Yf
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
        const label speciesIndex = thermo().species()[speciesName];

        return
            thermo().Y()[speciesIndex]
           *(scalar(1) - wRatioByP()*saturationModel_->pSat(Tf))
           /max(scalar(1) - thermo().Y()[saturatedIndex_], small);
    }
}


Foam::tmp<Foam::volScalarField>
Foam::interfaceCompositionModels::saturated::YfPrime
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
        const label speciesIndex = thermo().species()[speciesName];

        return
          - thermo().Y()[speciesIndex]
           *wRatioByP()*saturationModel_->pSatPrime(Tf)
           /max(scalar(1) - thermo().Y()[saturatedIndex_], small);
    }
}


// ************************************************************************* //
