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

SourceFiles
    tableReaders.C

\*---------------------------------------------------------------------------*/

#ifndef tableReaders_H
#define tableReaders_H

#include "tableReader.H"
#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Only used internally
#define makeTypeTableReadersTypeName(typeTableReader, dataType)               \
                                                                              \
    defineNamedTemplateTypeNameAndDebug(typeTableReader< dataType >, 0)

// Sometimes used externally
#define makeTableReadersTypeName(typeTableReader)                             \
                                                                              \
    makeTypeTableReadersTypeName(typeTableReader, scalar);                    \
    makeTypeTableReadersTypeName(typeTableReader, vector);                    \
    makeTypeTableReadersTypeName(typeTableReader, sphericalTensor);           \
    makeTypeTableReadersTypeName(typeTableReader, symmTensor);                \
    makeTypeTableReadersTypeName(typeTableReader, tensor)

// Define type info for single dataType template instantiation (eg, vector)
#define makeTableReaderType(typeTableReader, dataType)                        \
                                                                              \
    defineNamedTemplateTypeNameAndDebug(typeTableReader< dataType >, 0);      \
    addTemplatedToRunTimeSelectionTable                                       \
    (                                                                         \
        tableReader, typeTableReader, dataType, dictionary                    \
    )


// Define type info for scalar, vector etc. instantiations
#define makeTableReaders(typeTableReader)                                     \
                                                                              \
    makeTableReaderType(typeTableReader, scalar);                             \
    makeTableReaderType(typeTableReader, vector);                             \
    makeTableReaderType(typeTableReader, sphericalTensor);                    \
    makeTableReaderType(typeTableReader, symmTensor);                         \
    makeTableReaderType(typeTableReader, tensor)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
