/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "multiValveEngine.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshMovers::multiValveEngine::valveList::valveList
(
    const multiValveEngine& engine,
    const dictionary& dict
)
{
    // Count number of valves
    label nValves = 0;
    forAllConstIter(dictionary, dict, iter)
    {
        if (iter().isDict())
        {
            nValves++;
        }
    }

    // Add the valves to the list
    if (nValves > 0)
    {
        PtrList<valveObject>::setSize(nValves);

        label i = 0;
        forAllConstIter(dictionary, dict, iter)
        {
            if (iter().isDict())
            {
                const word& name = iter().keyword();
                const dictionary& valveDict = iter().dict();

                PtrList<valveObject>::set
                (
                    i++,
                    new valveObject(name, engine, valveDict)
                );
            }
        }
    }
}


// ************************************************************************* //
