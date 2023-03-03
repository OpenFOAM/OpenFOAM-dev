/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "fixedInterfacialAreaDiameter.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
    defineTypeNameAndDebug(fixedInterfacialArea, 0);
    addToRunTimeSelectionTable(diameterModel, fixedInterfacialArea, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::fixedInterfacialArea::fixedInterfacialArea
(
    const dictionary& diameterProperties,
    const phaseModel& phase
)
:
    diameterModel(diameterProperties, phase),
    AvbyAlpha_(inv(dimLength), -1),
    AvbyAlphaFieldPtr_()
{
    read(diameterProperties);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModels::fixedInterfacialArea::~fixedInterfacialArea()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::diameterModels::fixedInterfacialArea::d()
const
{
    if (AvbyAlphaFieldPtr_.valid())
    {
        const volScalarField& AvbyAlpha = AvbyAlphaFieldPtr_;
        return 6/AvbyAlpha;
    }

    return volScalarField::New
    (
        IOobject::groupName("d", phase().name()),
        phase().mesh(),
        6/AvbyAlpha_
    );
}


Foam::tmp<Foam::volScalarField> Foam::diameterModels::fixedInterfacialArea::Av()
const
{
    if (AvbyAlphaFieldPtr_.valid())
    {
        const volScalarField& AvbyAlpha = AvbyAlphaFieldPtr_;
        return phase()*AvbyAlpha;
    }

    return phase()*AvbyAlpha_;
}


bool Foam::diameterModels::fixedInterfacialArea::
read(const dictionary& phaseProperties)
{
    diameterModel::read(phaseProperties);

    AvbyAlpha_ = dimensionedScalar
        (
            inv(dimLength),
            diameterProperties().lookupOrDefault<scalar>("AvbyAlpha", -1)
        );

    if (AvbyAlpha_.value() < 0 && !AvbyAlphaFieldPtr_.valid())
    {
        Info<< "fixedInterfacialArea: Uniform AvbyAlpha not provided. "
            << "Looking up " << IOobject::groupName("AvbyAlpha", phase().name())
            << endl;

        AvbyAlphaFieldPtr_ =
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "AvbyAlpha",
                        phase().name()
                    ),
                    phase().mesh().time().name(),
                    phase().mesh(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                phase().mesh()
            );

        AvbyAlphaFieldPtr_->max
            (
                phaseProperties.lookupOrDefault<scalar>
                (
                    "residualAvbyAlpha",
                    rootSmall
                )
            );
    }

    return true;
}


// ************************************************************************* //
