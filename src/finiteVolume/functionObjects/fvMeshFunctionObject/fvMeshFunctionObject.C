/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2024 OpenFOAM Foundation
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

#include "fvMeshFunctionObject.H"
#include "Time.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fvMeshFunctionObject, 0);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

Foam::labelHashSet Foam::functionObjects::fvMeshFunctionObject::patchSet
(
    const dictionary& dict,
    const bool optional
) const
{
    if (dict.found("patch"))
    {
        const word patchName(dict.lookup("patch"));
        const label patchIndex = mesh_.boundaryMesh().findIndex(patchName);

        if (patchIndex >= 0)
        {
            return labelHashSet(FixedList<label, 1>(patchIndex));
        }
        else
        {
            FatalIOErrorInFunction(dict)
                << "Unable to find patch " << patchName << exit(FatalIOError);
        }
    }
    else if (dict.found("patches"))
    {
        return mesh_.boundaryMesh().patchSet
        (
            dict.lookup<wordReList>("patches")
        );
    }
    else
    {
        if (!optional)
        {
            FatalIOErrorInFunction(dict)
                << "Neither 'patch' or 'patches' specified"
                << exit(FatalIOError);
        }
    }

    return labelHashSet();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fvMeshFunctionObject::fvMeshFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    regionFunctionObject(name, runTime, dict),
    mesh_(refCast<const fvMesh>(obr_))
{}


Foam::functionObjects::fvMeshFunctionObject::fvMeshFunctionObject
(
    const word& name,
    const objectRegistry& obr
)
:
    regionFunctionObject(name, obr),
    mesh_(refCast<const fvMesh>(obr_))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fvMeshFunctionObject::~fvMeshFunctionObject()
{}


// ************************************************************************* //
