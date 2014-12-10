/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

Global
    Foam::ReadFields

Description
    Helper routine to read fields

SourceFiles
    ReadFields.C

\*---------------------------------------------------------------------------*/

#ifndef ReadFields_H
#define ReadFields_H

#include "PtrList.H"
#include "wordList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class IOobjectList;

//- Helper routine to read fields
template<class GeoField, class Mesh>
wordList ReadFields
(
    const Mesh& mesh,
    const IOobjectList& objects,
    PtrList<GeoField>& fields,
    const bool syncPar = true
);

//- Helper routine to read GeometricFields. The fieldsCache is per time
//  an objectRegistry of all stored fields
template<class GeoField>
static void ReadFields
(
    const word& fieldName,
    const typename GeoField::Mesh& mesh,
    const wordList& timeNames,
    objectRegistry& fieldsCache
);

//- Helper routine to read GeometricFields. The fieldsCache is per time
//  an objectRegistry of all stored fields
template<class GeoField>
static void ReadFields
(
    const word& fieldName,
    const typename GeoField::Mesh& mesh,
    const wordList& timeNames,
    const word& registryName = "fieldsCache"
);

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ReadFields.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
