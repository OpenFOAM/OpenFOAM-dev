/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::TableReader

Description
    Macros for creating standard TableReader-s

\*---------------------------------------------------------------------------*/

#ifndef makeTableReaders_H
#define makeTableReaders_H

#include "TableReader.H"
#include "EmbeddedTableReader.H"
#include "FoamTableReader.H"
#include "CsvTableReader.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeTableReaders(Type, nullArg)                                        \
                                                                               \
    defineTableReader(Type);                                                   \
                                                                               \
    namespace TableReaders                                                     \
    {                                                                          \
        addTableReader(Embedded, Type);                                        \
        addTableReader(Foam, Type);                                            \
        addTableReader(Csv, Type);                                             \
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
