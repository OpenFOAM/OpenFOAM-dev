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

#include "nonRandomTwoLiquid.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace interfaceCompositionModels
{
    defineTypeNameAndDebug(nonRandomTwoLiquid, 0);
    addToRunTimeSelectionTable
    (
        interfaceCompositionModel,
        nonRandomTwoLiquid,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceCompositionModels::nonRandomTwoLiquid::nonRandomTwoLiquid
(
    const dictionary& dict,
    const phaseInterface& interface
)
:
    interfaceCompositionModel(dict, interface),
    gamma1_
    (
        IOobject
        (
            IOobject::groupName("gamma1", this->interface().name()),
            interface.mesh().time().timeName(),
            interface.mesh()
        ),
        interface.mesh(),
        dimensionedScalar(dimless, 1)
    ),
    gamma2_
    (
        IOobject
        (
            IOobject::groupName("gamma2", this->interface().name()),
            interface.mesh().time().timeName(),
            interface.mesh()
        ),
        interface.mesh(),
        dimensionedScalar(dimless, 1)
    ),
    beta12_("", dimless/dimTemperature, 0),
    beta21_("", dimless/dimTemperature, 0)
{
    if (species().size() != 2)
    {
        FatalErrorInFunction
            << "nonRandomTwoLiquid model is suitable for two species only."
            << exit(FatalError);
    }

    species1Name_ = species()[0];
    species2Name_ = species()[1];

    species1Index_ = composition().species()[species1Name_];
    species2Index_ = composition().species()[species2Name_];

    alpha12_ = dimensionedScalar
    (
        "alpha12",
        dimless,
        dict.subDict(species1Name_).lookup("alpha")
    );
    alpha21_ = dimensionedScalar
    (
        "alpha21",
        dimless,
        dict.subDict(species2Name_).lookup("alpha")
    );

    beta12_ = dimensionedScalar
    (
        "beta12",
        dimless/dimTemperature,
        dict.subDict(species1Name_).lookup("beta")
    );
    beta21_ = dimensionedScalar
    (
        "beta21",
        dimless/dimTemperature,
        dict.subDict(species2Name_).lookup("beta")
    );

    saturationModel12_.reset
    (
        saturationModel::New
        (
            dict.subDict(species1Name_).subDict("interaction"),
            interface,
            false
        ).ptr()
    );
    saturationModel21_.reset
    (
        saturationModel::New
        (
            dict.subDict(species2Name_).subDict("interaction"),
            interface,
            false
        ).ptr()
    );

    speciesModel1_.reset
    (
        interfaceCompositionModel::New
        (
            dict.subDict(species1Name_),
            interface,
            false
        ).ptr()
    );
    speciesModel2_.reset
    (
        interfaceCompositionModel::New
        (
            dict.subDict(species2Name_),
            interface,
            false
        ).ptr()
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interfaceCompositionModels::nonRandomTwoLiquid::~nonRandomTwoLiquid()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::interfaceCompositionModels::nonRandomTwoLiquid::update
(
    const volScalarField& Tf
)
{
    const volScalarField W(thermo().W());

    const volScalarField X1
    (
        composition().Y(species1Index_)
       *W
       /dimensionedScalar
        (
            "W",
            dimMass/dimMoles,
            composition().Wi(species1Index_)
        )
    );

    const volScalarField X2
    (
        composition().Y(species2Index_)
       *W
       /dimensionedScalar
        (
            "W",
            dimMass/dimMoles,
            composition().Wi(species2Index_)
        )
    );

    const volScalarField alpha12(alpha12_ + Tf*beta12_);
    const volScalarField alpha21(alpha21_ + Tf*beta21_);

    const volScalarField tau12(saturationModel12_->lnPSat(Tf));
    const volScalarField tau21(saturationModel21_->lnPSat(Tf));

    const volScalarField G12(exp(- alpha12*tau12));
    const volScalarField G21(exp(- alpha21*tau21));

    gamma1_ =
        exp
        (
            sqr(X2)
           *(
                tau21*sqr(G21)/max(sqr(X1 + X2*G21), small)
              + tau12*G12/max(sqr(X2 + X1*G12), small)
            )
        );
    gamma2_ =
        exp
        (
            sqr(X1)
           *(
                tau12*sqr(G12)/max(sqr(X2 + X1*G12), small)
              + tau21*G21/max(sqr(X1 + X2*G21), small)
            )
        );
}


Foam::tmp<Foam::volScalarField>
Foam::interfaceCompositionModels::nonRandomTwoLiquid::Yf
(
    const word& speciesName,
    const volScalarField& Tf
) const
{
    if (speciesName == species1Name_)
    {
        return
            otherComposition().Y(speciesName)
           *speciesModel1_->Yf(speciesName, Tf)
           *gamma1_;
    }
    else if (speciesName == species2Name_)
    {
        return
            otherComposition().Y(speciesName)
           *speciesModel2_->Yf(speciesName, Tf)
           *gamma2_;
    }
    else
    {
        return
            composition().Y(speciesName)
           *(scalar(1) - Yf(species1Name_, Tf) - Yf(species2Name_, Tf));
    }
}


Foam::tmp<Foam::volScalarField>
Foam::interfaceCompositionModels::nonRandomTwoLiquid::YfPrime
(
    const word& speciesName,
    const volScalarField& Tf
) const
{
    if (speciesName == species1Name_)
    {
        return
            otherComposition().Y(speciesName)
           *speciesModel1_->YfPrime(speciesName, Tf)
           *gamma1_;
    }
    else if (speciesName == species2Name_)
    {
        return
            otherComposition().Y(speciesName)
           *speciesModel2_->YfPrime(speciesName, Tf)
           *gamma2_;
    }
    else
    {
        return
          - composition().Y(speciesName)
           *(YfPrime(species1Name_, Tf) + YfPrime(species2Name_, Tf));
    }
}


// ************************************************************************* //
