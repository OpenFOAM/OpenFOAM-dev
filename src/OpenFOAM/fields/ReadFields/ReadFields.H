/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Field reading functions for post-processing utilities

SourceFiles
    ReadFields.C

\*---------------------------------------------------------------------------*/

#ifndef ReadFields_H
#define ReadFields_H

#include "PtrList.H"
#include "wordList.H"
#include "HashSet.H"
#include "LIFOStack.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class regIOobject;
class IOobjectList;
class objectRegistry;

//- Read all fields of the specified type.
//  Returns names of fields read.
//  Guarantees all processors read fields in same order.
template<class GeoField, class Mesh>
wordList ReadFields
(
    const Mesh& mesh,
    const IOobjectList& objects,
    PtrList<GeoField>& fields,
    const bool syncPar = true
);

//- Read all GeometricFields of the specified type.
//  The fieldsCache is an objectRegistry of all stored fields
template<class GeoField>
static void ReadFields
(
    const word& fieldName,
    const typename GeoField::Mesh& mesh,
    const wordList& timeNames,
    objectRegistry& fieldsCache
);

//- Read all GeometricFields of the specified type.
//  The fieldsCache is an objectRegistry of all stored fields
template<class GeoField>
static void ReadFields
(
    const word& fieldName,
    const typename GeoField::Mesh& mesh,
    const wordList& timeNames,
    const word& registryName = "fieldsCache"
);

//- Read the selected GeometricFields of the specified type.
//  The fields are transferred to the objectRegistry and a list of them is
//  returned as a stack for later clean-up
template<class GeoFieldType>
void readFields
(
    const typename GeoFieldType::Mesh& mesh,
    const IOobjectList& objects,
    const HashSet<word>& selectedFields,
    LIFOStack<regIOobject*>& storedObjects
);


//- Read the selected UniformDimensionedFields of the specified type.
//  The fields are transferred to the objectRegistry and a list of them is
//  returned as a stack for later clean-up
template<class GeoFieldType>
void readUniformFields
(
    const IOobjectList& objects,
    const HashSet<word>& selectedFields,
    LIFOStack<regIOobject*>& storedObjects,
    const bool syncPar = true
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ReadFields.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
