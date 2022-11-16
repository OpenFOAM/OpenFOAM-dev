/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2022 OpenFOAM Foundation
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

#include "cavitation.H"
#include "phaseSystem.H"
#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseTransferModels
{
    defineTypeNameAndDebug(cavitation, 0);
    addToRunTimeSelectionTable
    (
        phaseTransferModel,
        cavitation,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseTransferModels::cavitation::cavitation
(
    const dictionary& dict,
    const phaseInterface& interface
)
:
    phaseTransferModel(dict, interface),
    interface_(interface),
    cavitation_(compressible::cavitationModel::New(dict, interface_))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseTransferModels::cavitation::~cavitation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::phaseTransferModels::cavitation::mixture() const
{
    return true;
}


Foam::tmp<Foam::volScalarField>
Foam::phaseTransferModels::cavitation::dmdtf() const
{
    tmp<volScalarField> tResult =
        volScalarField::New
        (
            IOobject::groupName(typedName("dmdtf"), interface_.name()),
            interface_.mesh(),
            dimDensity/dimTime,
            zeroGradientFvPatchField<scalar>::typeName
        );

    const Pair<tmp<volScalarField::Internal>> coeffs(cavitation_->mDot12P());

    const volScalarField::Internal& p = interface_.phase1().thermo().p();
    const volScalarField::Internal pSat1(cavitation_->pSat1());
    const volScalarField::Internal pSat2(cavitation_->pSat2());

    tResult.ref().ref() = coeffs[0]*(p - pSat1) - coeffs[1]*(p - pSat2);
    tResult.ref().correctBoundaryConditions();

    return tResult;
}


Foam::tmp<Foam::volScalarField>
Foam::phaseTransferModels::cavitation::d2mdtdpf() const
{
    tmp<volScalarField> tResult =
        volScalarField::New
        (
            IOobject::groupName(typedName("d2mdtdpf"), interface_.name()),
            interface_.mesh(),
            dimDensity/dimTime/dimPressure,
            zeroGradientFvPatchField<scalar>::typeName
        );

    const Pair<tmp<volScalarField::Internal>> coeffs(cavitation_->mDot12P());

    tResult.ref().ref() = coeffs[0] - coeffs[1];
    tResult.ref().correctBoundaryConditions();

    return tResult;
}


// ************************************************************************* //
