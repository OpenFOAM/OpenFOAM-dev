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

#include "phaseSolidThermophysicalTransportModel.H"
#include "isotropic.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineRunTimeSelectionTable
    (
        phaseSolidThermophysicalTransportModel,
        dictionary
    );
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::phaseSolidThermophysicalTransportModel::printCoeffs
(const word& type)
{
    if (printCoeffs_)
    {
        Info<< coeffDict_.dictName() << coeffDict_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseSolidThermophysicalTransportModel::
phaseSolidThermophysicalTransportModel
(
    const word& type,
    const alphaField& alpha,
    const solidThermo& thermo
)
:
    thermophysicalTransportModel(thermo.mesh(), alpha.group()),
    alpha_(alpha),
    thermo_(thermo),
    printCoeffs_(lookupOrDefault<Switch>("printCoeffs", false)),
    coeffDict_(optionalSubDict(type + "Coeffs"))
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::phaseSolidThermophysicalTransportModel>
Foam::phaseSolidThermophysicalTransportModel::New
(
    const alphaField& alpha,
    const solidThermo& thermo
)
{
    typeIOobject<IOdictionary> header
    (
        IOobject::groupName
        (
            phaseSolidThermophysicalTransportModel::typeName,
            alpha.group()
        ),
        thermo.mesh().time().constant(),
        thermo.mesh(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    if (header.headerOk())
    {
        const word modelType(IOdictionary(header).lookup("model"));

        Info<< "Selecting solid thermophysical transport model "
            << modelType << endl;

        typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            FatalErrorInFunction
                << "Unknown solid thermophysical transport model "
                << modelType << nl << nl
                << "Available models:" << endl
                << dictionaryConstructorTablePtr_->sortedToc()
                << exit(FatalError);
        }

        return autoPtr<phaseSolidThermophysicalTransportModel>
        (
            cstrIter()(alpha, thermo)
        );
    }
    else
    {
        Info<< "Selecting default solid thermophysical transport model "
            << solidThermophysicalTransportModels::
               isotropic<phaseSolidThermophysicalTransportModel>::typeName
            << endl;

        return autoPtr<phaseSolidThermophysicalTransportModel>
        (
            new solidThermophysicalTransportModels::
            isotropic<phaseSolidThermophysicalTransportModel>(alpha, thermo)
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::phaseSolidThermophysicalTransportModel::kappa() const
{
    return thermo().kappa();
}


Foam::tmp<Foam::scalarField>
Foam::phaseSolidThermophysicalTransportModel::kappa
(
    const label patchi
) const
{
    return thermo().kappa().boundaryField()[patchi];
}


bool Foam::phaseSolidThermophysicalTransportModel::read()
{
    if (regIOobject::read())
    {
        coeffDict_ <<= optionalSubDict(type() + "Coeffs");
        return true;
    }
    else
    {
        return false;
    }
}


void Foam::phaseSolidThermophysicalTransportModel::predict()
{}


void Foam::phaseSolidThermophysicalTransportModel::correct()
{}


// ************************************************************************* //
