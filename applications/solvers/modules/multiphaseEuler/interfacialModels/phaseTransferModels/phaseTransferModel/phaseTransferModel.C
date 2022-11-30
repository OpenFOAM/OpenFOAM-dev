/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2022 OpenFOAM Foundation
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

#include "phaseTransferModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseTransferModel, 0);
    defineBlendedInterfacialModelTypeNameAndDebug(phaseTransferModel, 0);
    defineRunTimeSelectionTable(phaseTransferModel, dictionary);
}

const Foam::dimensionSet Foam::phaseTransferModel::dimDmdt =
    Foam::dimDensity/Foam::dimTime;

const Foam::dimensionSet Foam::phaseTransferModel::dimD2mdtdp =
    Foam::dimDensity/Foam::dimTime/Foam::dimPressure;

const Foam::hashedWordList Foam::phaseTransferModel::noSpecies_ =
    Foam::hashedWordList();

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseTransferModel::phaseTransferModel
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
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseTransferModel::~phaseTransferModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::phaseTransferModel::mixture() const
{
    return false;
}


Foam::tmp<Foam::volScalarField> Foam::phaseTransferModel::dmdtf() const
{
    return tmp<volScalarField>(nullptr);
}


Foam::tmp<Foam::volScalarField> Foam::phaseTransferModel::d2mdtdpf() const
{
    return tmp<volScalarField>(nullptr);
}


const Foam::hashedWordList& Foam::phaseTransferModel::species() const
{
    return noSpecies_;
}


Foam::HashPtrTable<Foam::volScalarField>
Foam::phaseTransferModel::dmidtf() const
{
    return HashPtrTable<volScalarField>();
}


bool Foam::phaseTransferModel::writeData(Ostream& os) const
{
    return os.good();
}


bool Foam::blendedPhaseTransferModel::mixture() const
{
    return evaluate(&phaseTransferModel::mixture);
}


Foam::tmp<Foam::volScalarField>
Foam::blendedPhaseTransferModel::dmdtf() const
{
    return
        evaluate
        (
            &phaseTransferModel::dmdtf,
            "dmdtf",
            phaseTransferModel::dimDmdt,
            true
        );
}


Foam::tmp<Foam::volScalarField>
Foam::blendedPhaseTransferModel::d2mdtdpf() const
{
    return
        evaluate
        (
            &phaseTransferModel::d2mdtdpf,
            "d2mdtdpf",
            phaseTransferModel::dimD2mdtdp,
            true
        );
}


Foam::hashedWordList Foam::blendedPhaseTransferModel::species() const
{
    return evaluate(&phaseTransferModel::species);
}


Foam::HashPtrTable<Foam::volScalarField>
Foam::blendedPhaseTransferModel::dmidtf() const
{
    return
        evaluate
        (
            &phaseTransferModel::dmidtf,
            "dmidtf",
            phaseTransferModel::dimDmdt,
            true
        );
}


// ************************************************************************* //
