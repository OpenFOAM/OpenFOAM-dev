/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "coordinateSystems.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace coordinateSystems
{
    defineTypeNameAndDebug(coordinateSystems, 0);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coordinateSystems::coordinateSystems::coordinateSystems
(
    const objectRegistry& obr
)
:
    MeshObject<objectRegistry, GeometricMeshObject, coordinateSystems>
    (
        obr,
        IOobject
        (
            typeName,
            obr.time().constant(),
            obr,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    PtrDictionary<coordinateSystem>()
{
    readHeaderOk(IOstream::ASCII, typeName);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::coordinateSystems::coordinateSystems::readData(Istream& is)
{
    const dictionary coordinateSystemsDict(is);

    forAllConstIter(dictionary, coordinateSystemsDict, iter)
    {
        if (iter().isDict())
        {
            const word& name = iter().keyword();
            const dictionary& dict = iter().dict();

            this->insert
            (
                name,
                coordinateSystem::New(name, dict).ptr()
            );
        }
    }

    return !is.bad();
}


// ************************************************************************* //
