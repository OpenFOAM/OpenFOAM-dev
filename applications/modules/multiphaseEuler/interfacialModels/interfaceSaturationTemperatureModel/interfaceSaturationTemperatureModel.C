/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2025 OpenFOAM Foundation
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

#include "interfaceSaturationTemperatureModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(interfaceSaturationTemperatureModel, 0);
    defineRunTimeSelectionTable
    (
        interfaceSaturationTemperatureModel,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceSaturationTemperatureModel::interfaceSaturationTemperatureModel
(
    const dictionary& dict,
    const phaseInterface& interface
)
:
    regIOobject
    (
        IOobject
        (
            IOobject::groupName(typeName, interface.name()),
            interface.mesh().time().name(),
            interface.mesh()
        )
    ),
    saturationModel_(saturationTemperatureModel::New(dict)),
    interface_(interface)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interfaceSaturationTemperatureModel::
~interfaceSaturationTemperatureModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::phaseInterface&
Foam::interfaceSaturationTemperatureModel::interface() const
{
    return interface_;
}


Foam::tmp<Foam::scalarField>
Foam::interfaceSaturationTemperatureModel::Tsat
(
    const scalarField& p
) const
{
    return saturationModel_->Tsat(p);
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::interfaceSaturationTemperatureModel::Tsat
(
    const volScalarField::Internal& p
) const
{
    return saturationModel_->Tsat(p);
}


Foam::tmp<Foam::volScalarField>
Foam::interfaceSaturationTemperatureModel::Tsat
(
    const volScalarField& p
) const
{
    return saturationModel_->Tsat(p);
}


Foam::tmp<Foam::scalarField>
Foam::interfaceSaturationTemperatureModel::TsatPrime
(
    const scalarField& p
) const
{
    return saturationModel_->TsatPrime(p);
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::interfaceSaturationTemperatureModel::TsatPrime
(
    const volScalarField::Internal& p
) const
{
    return saturationModel_->TsatPrime(p);
}


Foam::tmp<Foam::volScalarField>
Foam::interfaceSaturationTemperatureModel::TsatPrime
(
    const volScalarField& p
) const
{
    return saturationModel_->TsatPrime(p);
}


bool Foam::interfaceSaturationTemperatureModel::writeData(Ostream& os) const
{
    return os.good();
}


// ************************************************************************* //
