/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "polyTopoChanger.H"
#include "polyMesh.H"
#include "polyTopoChange.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(polyTopoChanger, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polyTopoChanger::polyTopoChanger(polyMesh& mesh)
:
    PtrList<polyMeshModifier>(),
    mesh_(mesh)
{}


bool Foam::polyTopoChanger::changeTopology() const
{
    // Go through all mesh modifiers and accumulate the morphing information
    const PtrList<polyMeshModifier>& topoChanges = *this;

    bool triggerChange = false;

    forAll(topoChanges, morphI)
    {
        bool curTriggerChange = topoChanges[morphI].changeTopology();

        if (debug)
        {
            Info<< "Modifier " << morphI << " named "
                << topoChanges[morphI].name();

            if (curTriggerChange)
            {
                Info<< " morphing" << endl;
            }
            else
            {
                Info<< " unchanged" << endl;
            }
        }

        triggerChange = triggerChange || curTriggerChange;
    }

    return triggerChange;
}


Foam::autoPtr<Foam::polyTopoChange>
Foam::polyTopoChanger::topoChangeRequest() const
{
    // Collect changes from all modifiers
    const PtrList<polyMeshModifier>& topoChanges = *this;

    polyTopoChange* refPtr(new polyTopoChange(mesh()));
    polyTopoChange& ref = *refPtr;

    forAll(topoChanges, morphI)
    {
        topoChanges[morphI].setRefinement(ref);
    }

    return autoPtr<polyTopoChange>(refPtr);
}


void Foam::polyTopoChanger::update(const polyTopoChangeMap& m)
{
    // Go through all mesh modifiers and accumulate the morphing information
    PtrList<polyMeshModifier>& topoChanges = *this;

    forAll(topoChanges, morphI)
    {
        topoChanges[morphI].topoChange(m);
    }
}


Foam::autoPtr<Foam::polyTopoChangeMap> Foam::polyTopoChanger::changeMesh
(
    const bool inflate,
    const bool syncParallel,
    const bool orderCells,
    const bool orderPoints
)
{
    if (changeTopology())
    {
        autoPtr<polyTopoChange> ref = topoChangeRequest();

        autoPtr<polyTopoChangeMap> topoChangeMap = ref().changeMesh
        (
            mesh_,
            inflate,
            syncParallel,
            orderCells,
            orderPoints
        );

        update(topoChangeMap());
        mesh_.topoChange(topoChangeMap());
        return topoChangeMap;
    }
    else
    {
        return autoPtr<polyTopoChangeMap>(nullptr);
    }
}


// ************************************************************************* //
