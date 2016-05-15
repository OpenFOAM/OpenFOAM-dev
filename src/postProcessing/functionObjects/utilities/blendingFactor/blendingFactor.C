/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
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

#include "blendingFactor.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(blendingFactor, 0);
    addToRunTimeSelectionTable(functionObject, blendingFactor, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::blendingFactor::blendingFactor
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObject(name),
    obr_
    (
        runTime.lookupObject<objectRegistry>
        (
            dict.lookupOrDefault("region", polyMesh::defaultRegion)
        )
    ),
    phiName_("unknown-phiName"),
    fieldName_("unknown-fieldName")
{
    if (!isA<fvMesh>(obr_))
    {
        FatalErrorInFunction
            << "objectRegistry is not an fvMesh" << exit(FatalError);
    }

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::blendingFactor::~blendingFactor()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::blendingFactor::read(const dictionary& dict)
{
    phiName_ = dict.lookupOrDefault<word>("phiName", "phi");
    dict.lookup("fieldName") >> fieldName_;

    return true;
}


bool Foam::functionObjects::blendingFactor::execute(const bool postProcess)
{
    calc<scalar>();
    calc<vector>();

    return true;
}


bool Foam::functionObjects::blendingFactor::write(const bool postProcess)
{
    const word fieldName = "blendingFactor:" + fieldName_;

    const volScalarField& blendingFactor =
        obr_.lookupObject<volScalarField>(fieldName);

    Info<< type() << " " << name() << " output:" << nl
        << "    writing field " << blendingFactor.name() << nl
        << endl;

    blendingFactor.write();

    return true;
}


// ************************************************************************* //
