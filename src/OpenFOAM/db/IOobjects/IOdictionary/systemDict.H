/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2024 OpenFOAM Foundation
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
    Foam::systemDict

Description
    Read and return the specified dictionary from system or from path provided
    with the -dict option.

SourceFiles
    systemDict.C

\*---------------------------------------------------------------------------*/

#ifndef systemDict_H
#define systemDict_H

#include "argList.H"
#include "IOdictionary.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * * //

IOobject systemDictIO
(
    const word& dictName,
    const argList& args,
    const objectRegistry& ob,
    const word& regionName = polyMesh::defaultRegion,
    const fileName& path = fileName::null
);

IOdictionary systemDict
(
    const word& dictName,
    const argList& args,
    const objectRegistry& ob,
    const word& regionName = polyMesh::defaultRegion,
    const fileName& path = fileName::null
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
