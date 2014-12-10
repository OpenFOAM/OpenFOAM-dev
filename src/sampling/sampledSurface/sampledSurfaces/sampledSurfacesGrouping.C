/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "sampledSurfaces.H"
#include "volFields.H"
#include "IOobjectList.H"
#include "stringListOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::sampledSurfaces::classifyFields()
{
     // check files for a particular time
    if (loadFromFiles_)
    {
        IOobjectList objects(mesh_, mesh_.time().timeName());
        wordList allFields = objects.sortedNames();
        labelList indices = findStrings(fieldSelection_, allFields);
        return indices.size();
    }
    else
    {
        wordList allFields = mesh_.sortedNames();
        labelList indices = findStrings(fieldSelection_, allFields);
        return indices.size();
    }
}


// ************************************************************************* //
