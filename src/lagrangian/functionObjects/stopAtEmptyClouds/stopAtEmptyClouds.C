/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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

#include "stopAtEmptyClouds.H"
#include "parcelCloudList.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(stopAtEmptyClouds, 0);
    addToRunTimeSelectionTable(functionObject, stopAtEmptyClouds, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::stopAtEmptyClouds::condition() const
{
    // Potentially nothing has been injected yet. Wait for a step. If the
    // clouds begin injection later than this, then time controls will be
    // needed to delay the execution of this function.
    if (time_.timeIndex() == time_.startTimeIndex())
    {
        return false;
    }

    // Loop all regions and all clouds in each region. If a non-empty cloud is
    // found then do not stop.
    bool noClouds = true;

    const HashTable<const polyMesh*> meshes =
        time_.lookupClass<polyMesh>();

    forAllConstIter(HashTable<const polyMesh*>, meshes, meshIter)
    {
        const HashTable<const parcelCloud*> clouds =
            meshIter()->lookupClass<parcelCloud>();

        forAllConstIter(HashTable<const parcelCloud*>, clouds, cloudIter)
        {
            noClouds = false;

            if (returnReduce(cloudIter()->nParcels(), sumOp<label>()) != 0)
            {
                return false;
            }
        }
    }

    // Either there are no clouds or there are clouds and all of them are
    // empty. Stop if the latter.
    return !noClouds;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::stopAtEmptyClouds::stopAtEmptyClouds
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    stopAt(name, runTime, dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::stopAtEmptyClouds::~stopAtEmptyClouds()
{}


// ************************************************************************* //
