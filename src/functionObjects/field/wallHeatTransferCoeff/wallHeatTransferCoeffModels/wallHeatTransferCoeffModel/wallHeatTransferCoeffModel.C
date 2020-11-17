/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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

#include "wallHeatTransferCoeffModel.H"
#include "fluidThermoMomentumTransportModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(wallHeatTransferCoeffModel, 0);
    defineRunTimeSelectionTable(wallHeatTransferCoeffModel, word);
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::wallHeatTransferCoeffModel>
Foam::wallHeatTransferCoeffModel::New
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
{
    const word model(dict.lookupOrDefault<word>("model", "kappaEff"));

    if (debug)
    {
        Info<< "Selecting heat transfer coefficient type: "
            << model << endl;
    }

    wordConstructorTable::iterator cstrIter =
        wordConstructorTablePtr_->find(model);

    if (cstrIter == wordConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown heat transfer coefficient type "
            << model << nl << nl
            << "Valid coefficient types: " << endl
            << wordConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<wallHeatTransferCoeffModel>(cstrIter()(name, mesh, dict));
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volSymmTensorField> Foam::wallHeatTransferCoeffModel::tau
(
    const momentumTransportModel& mmtm,
    const fvMesh& mesh
) const
{
    if (isA<incompressible::momentumTransportModel>(mmtm))
    {
        return
            refCast<const incompressible::momentumTransportModel>
            (
                mmtm
            ).devSigma();
    }
    else if (isA<compressible::momentumTransportModel>(mmtm))
    {
        return
            refCast<const compressible::momentumTransportModel>(mmtm).devTau()
           /refCast<const compressible::momentumTransportModel>(mmtm).rho();
    }
    else
    {
        FatalErrorInFunction
            << "The type of momentum transport model was not recognised"
            << exit(FatalError);
    }

    return tmp<Foam::volSymmTensorField>();
}

// ************************************************************************* //
