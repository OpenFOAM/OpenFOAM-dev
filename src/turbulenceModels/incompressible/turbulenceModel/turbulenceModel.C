/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "turbulenceModel.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(turbulenceModel, 0);
defineRunTimeSelectionTable(turbulenceModel, turbulenceModel);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulenceModel::turbulenceModel
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName
)
:
    regIOobject
    (
        IOobject
        (
            turbulenceModelName,
            U.time().constant(),
            U.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    runTime_(U.time()),
    mesh_(U.mesh()),

    U_(U),
    phi_(phi),
    transportModel_(transport),
    y_(mesh_)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<turbulenceModel> turbulenceModel::New
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName
)
{
    // get model name, but do not register the dictionary
    // otherwise it is registered in the database twice
    const word modelType
    (
        IOdictionary
        (
            IOobject
            (
                "turbulenceProperties",
                U.time().constant(),
                U.db(),
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        ).lookup("simulationType")
    );

    Info<< "Selecting turbulence model type " << modelType << endl;

    turbulenceModelConstructorTable::iterator cstrIter =
        turbulenceModelConstructorTablePtr_->find(modelType);

    if (cstrIter == turbulenceModelConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "turbulenceModel::New(const volVectorField&, "
            "const surfaceScalarField&, transportModel&, const word&)"
        )   << "Unknown turbulenceModel type "
            << modelType << nl << nl
            << "Valid turbulenceModel types:" << endl
            << turbulenceModelConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<turbulenceModel>
    (
        cstrIter()(U, phi, transport, turbulenceModelName)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void turbulenceModel::correct()
{
    transportModel_.correct();

    if (mesh_.changing())
    {
        y_.correct();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
