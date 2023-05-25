/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "heatTransferModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(heatTransferModel, 0);
    defineBlendedInterfacialModelTypeNameAndDebug(heatTransferModel, 0);
    defineSidedInterfacialModelTypeNameAndDebug(blendedHeatTransferModel, 0);
    defineRunTimeSelectionTable(heatTransferModel, dictionary);
}

const Foam::dimensionSet Foam::heatTransferModel::dimK(1, -1, -3, -1, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatTransferModel::heatTransferModel
(
    const dictionary& dict,
    const phaseInterface& interface,
    const bool registerObject
)
:
    regIOobject
    (
        IOobject
        (
            IOobject::groupName(typeName, interface.name()),
            interface.mesh().time().name(),
            interface.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            registerObject
        )
    ),
    residualAlpha_
    (
        "residualAlpha",
        dimless,
        dict.lookupOrDefault<scalar>
        (
            "residualAlpha",
            sqrt
            (
                interface.phase1().residualAlpha().value()
               *interface.phase2().residualAlpha().value()
            )
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::heatTransferModel::~heatTransferModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::heatTransferModel::K() const
{
    return K(residualAlpha_.value());
}


bool Foam::heatTransferModel::writeData(Ostream& os) const
{
    return os.good();
}


Foam::tmp<Foam::volScalarField> Foam::blendedHeatTransferModel::K() const
{
    tmp<volScalarField> (heatTransferModel::*k)() const =
        &heatTransferModel::K;
    return evaluate(k, "K", heatTransferModel::dimK, false);
}


Foam::tmp<Foam::volScalarField> Foam::blendedHeatTransferModel::K
(
    const scalar residualAlpha
) const
{
    tmp<volScalarField> (heatTransferModel::*k)(const scalar) const =
        &heatTransferModel::K;
    return evaluate(k, "Kf", heatTransferModel::dimK, false, residualAlpha);
}


// ************************************************************************* //
