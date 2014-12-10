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

InNamespace
    Foam::mapClouds

Description
    Generic Geometric field mapper.  For "real" mapping, add template
    specialisations for mapping of internal fields depending on mesh
    type.

\*---------------------------------------------------------------------------*/

#ifndef mapClouds_H
#define mapClouds_H

#include "cloud.H"
#include "objectRegistry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Generic Geometric field mapper.
//  For "real" mapping, add template specialisations for mapping of internal
//  fields depending on mesh type.
inline void mapClouds(const objectRegistry& db, const mapPolyMesh& mapper)
{
    HashTable<const cloud*> clouds(db.lookupClass<cloud>());

    forAllIter(HashTable<const cloud*>, clouds, iter)
    {
        cloud& c = const_cast<cloud&>(*iter());

        if (polyMesh::debug)
        {
            Info<< "Mapping cloud " << c.name() << endl;
        }

        c.autoMap(mapper);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
