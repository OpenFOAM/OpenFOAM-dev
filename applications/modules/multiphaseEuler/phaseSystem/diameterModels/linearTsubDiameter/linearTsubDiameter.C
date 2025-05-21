/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2025 OpenFOAM Foundation
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

#include "linearTsubDiameter.H"
#include "phaseSystem.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
    defineTypeNameAndDebug(linearTsub, 0);
    addToRunTimeSelectionTable(diameterModel, linearTsub, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::linearTsub::linearTsub
(
    const dictionary& diameterProperties,
    const phaseModel& phase
)
:
    spherical(diameterProperties, phase),
    liquidPhaseName_(diameterProperties.lookup("liquidPhase")),
    d2_("d2", dimLength, diameterProperties.lookupOrDefault("d2", 0.0015)),
    Tsub2_
    (
        "Tsub2",
         dimTemperature,
         diameterProperties.lookupOrDefault("Tsub2", 0)
    ),
    d1_("d1", dimLength, diameterProperties.lookupOrDefault("d1", 0.00015)),
    Tsub1_
    (
        "Tsub1",
        dimTemperature,
        diameterProperties.lookupOrDefault("Tsub1", 13.5)
    ),
    saturationModelPtr_
    (
        saturationTemperatureModel::New
        (
            "saturationTemperature",
            diameterProperties
        ).ptr()
    ),
    d_
    (
        IOobject
        (
            IOobject::groupName("d", phase.name()),
            phase.time().name(),
            phase.mesh()
        ),
        phase.mesh(),
        d1_
    )
{
    Info<< "    d2: " << d2_.value() << endl
        << "    Tsub2: " << Tsub2_.value() << endl
        << "    d1: " << d1_.value() << endl
        << "    Tsub1: " << Tsub1_.value() << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModels::linearTsub::~linearTsub()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::diameterModels::linearTsub::d() const
{
    return d_;
}


void Foam::diameterModels::linearTsub::correct()
{
    const phaseModel& liquid = phase().fluid().phases()[liquidPhaseName_];

    const volScalarField Tsub
    (
        saturationModelPtr_->Tsat(liquid.fluidThermo().p())
      - liquid.thermo().T()
    );

    d_ = max
    (
        d1_,
        min
        (
            d2_,
            (d1_*(Tsub - Tsub2_) + d2_*(Tsub - Tsub1_))/(Tsub2_ - Tsub1_)
        )
    );
}


bool Foam::diameterModels::linearTsub::read(const dictionary& phaseProperties)
{
    spherical::read(phaseProperties);

    diameterProperties().lookup("liquidPhase") >> liquidPhaseName_;

    d2_.readIfPresent(diameterProperties());
    Tsub2_.readIfPresent(diameterProperties());
    d1_.readIfPresent(diameterProperties());
    Tsub1_.readIfPresent(diameterProperties());

    saturationModelPtr_.reset
    (
        saturationTemperatureModel::New
        (
            "saturationTemperature",
            diameterProperties()
        ).ptr()
    );

    return true;
}


// ************************************************************************* //
