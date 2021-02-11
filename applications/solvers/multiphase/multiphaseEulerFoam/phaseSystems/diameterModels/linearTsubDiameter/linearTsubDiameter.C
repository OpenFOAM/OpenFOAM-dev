/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2021 OpenFOAM Foundation
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
#include "saturationModel.H"
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
    d1_
    (
        "d1",
        dimLength,
        diameterProperties.lookupOrDefault("d1", 0.00015)
    ),
    Tsub1_
    (
        "Tsub1",
        dimTemperature,
        diameterProperties.lookupOrDefault("Tsub1", 13.5)
    ),
    d_
    (
        IOobject
        (
            IOobject::groupName("d", phase.name()),
            phase.time().timeName(),
            phase.mesh()
        ),
        phase.mesh(),
        d1_
    )
{}


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
    // Lookup the fluid model
    const phaseSystem& fluid =
        refCast<const phaseSystem>
        (
            phase().mesh().lookupObject<phaseSystem>("phaseProperties")
        );

    const phaseModel& liquid(fluid.phases()[liquidPhaseName_]);

    if (phase().mesh().foundObject<saturationModel>("saturationModel"))
    {
        const saturationModel& satModel =
            phase().mesh().lookupObject<saturationModel>("saturationModel");

        const volScalarField Tsub
        (
            liquid.thermo().T() - satModel.Tsat(liquid.thermo().p())
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
}


bool Foam::diameterModels::linearTsub::read(const dictionary& phaseProperties)
{
    spherical::read(phaseProperties);

    diameterProperties().lookup("liquidPhase") >> liquidPhaseName_;
    diameterProperties().lookup("d2") >> d2_;
    diameterProperties().lookup("Tsub2") >> Tsub2_;
    diameterProperties().lookup("d1") >> d1_;
    diameterProperties().lookup("Tsub1") >> Tsub1_;

    return true;
}


// ************************************************************************* //
