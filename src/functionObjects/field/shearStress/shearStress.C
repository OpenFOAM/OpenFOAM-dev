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

#include "shearStress.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(shearStress, 0);
    addToRunTimeSelectionTable(functionObject, shearStress, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::shearStress::shearStress
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeLocalObjects(obr_, false),
    phaseName_(word::null)
{
    read(dict);
    resetLocalObjectName(IOobject::groupName(type(), phaseName_));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::shearStress::~shearStress()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::shearStress::read
(
    const dictionary& dict
)
{
    fvMeshFunctionObject::read(dict);
    writeLocalObjects::read(dict);

    phaseName_ = dict.lookupOrDefault<word>("phase", word::null);

    return true;
}


bool Foam::functionObjects::shearStress::execute()
{
    const word fieldName(IOobject::groupName(type(), phaseName_));

    typedef compressibleTurbulenceModel cmpModel;
    typedef incompressibleTurbulenceModel icoModel;

    const word turbulenceModelName
    (
        IOobject::groupName(turbulenceModel::propertiesName, phaseName_)
    );

    if (mesh_.foundObject<cmpModel>(turbulenceModelName))
    {
        const cmpModel& model =
            mesh_.lookupObject<cmpModel>(turbulenceModelName);

        return store(fieldName, model.devRhoReff());
    }
    else if (mesh_.foundObject<icoModel>(turbulenceModelName))
    {
        const icoModel& model =
            mesh_.lookupObject<icoModel>(turbulenceModelName);

        return store(fieldName, model.devReff());
    }
    else
    {
        FatalErrorInFunction
            << "Unable to find compressible turbulence model "
            << turbulenceModelName << " in the database"
            << exit(FatalError);

        return false;
    }
}


bool Foam::functionObjects::shearStress::write()
{
    return writeLocalObjects::write();
}


// ************************************************************************* //
