/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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

#include "physicalProperties.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(physicalProperties, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * //

Foam::IOobject Foam::physicalProperties::findModelDict
(
    const objectRegistry& obr,
    const word& group,
    bool registerObject
)
{
    typeIOobject<IOdictionary> physicalPropertiesIO
    (
        IOobject::groupName(physicalProperties::typeName, group),
        obr.time().constant(),
        obr,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE,
        registerObject
    );

    if (physicalPropertiesIO.headerOk())
    {
        return physicalPropertiesIO;
    }
    else
    {
        typeIOobject<IOdictionary> thermophysicalPropertiesIO
        (
            IOobject::groupName("thermophysicalProperties", group),
            obr.time().constant(),
            obr,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            registerObject
        );

        if (thermophysicalPropertiesIO.headerOk())
        {
            return thermophysicalPropertiesIO;
        }
        else
        {
            typeIOobject<IOdictionary> transportPropertiesIO
            (
                IOobject::groupName("transportProperties", group),
                obr.time().constant(),
                obr,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                registerObject
            );

            if (transportPropertiesIO.headerOk())
            {
                return transportPropertiesIO;
            }
            else
            {
                return physicalPropertiesIO;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::physicalProperties::physicalProperties
(
    const fvMesh& mesh,
    const word& group
)
:
    IOdictionary(findModelDict(mesh, group, true))
{
    // Ensure name of IOdictionary is typeName
    rename(IOobject::groupName(physicalProperties::typeName, group));
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::physicalProperties::read()
{
    return regIOobject::read();
}


// ************************************************************************* //
