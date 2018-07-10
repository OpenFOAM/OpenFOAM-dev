/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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

Class
    Foam::simpleObjectRegistry

Description
    Object registry for simpleRegIOobject. Maintains ordering.

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef simpleObjectRegistry_H
#define simpleObjectRegistry_H

#include "Dictionary.H"
#include "simpleRegIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class simpleObjectRegistryEntry Declaration
\*---------------------------------------------------------------------------*/

class simpleObjectRegistryEntry
:
    public Dictionary<simpleObjectRegistryEntry>::link,
    public List<simpleRegIOobject*>
{
public:

    simpleObjectRegistryEntry(const List<simpleRegIOobject*>& data)
    :
        List<simpleRegIOobject*>(data)
    {}
};


/*---------------------------------------------------------------------------*\
                      Class simpleObjectRegistry Declaration
\*---------------------------------------------------------------------------*/

class simpleObjectRegistry
:
    public Dictionary<simpleObjectRegistryEntry>
{
public:

    // Constructors

        //- Construct given initial table size
        simpleObjectRegistry(const label nIoObjects = 128)
        :
            Dictionary<simpleObjectRegistryEntry>(nIoObjects)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
