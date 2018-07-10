/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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

Description
    General C-preprocessor macros

\*---------------------------------------------------------------------------*/

#ifndef macros_H
#define macros_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Concatenate two preprocessor tokens
#define CAT_(a, b) a ## b
#define CAT(a, b) CAT_(a, b)

//- Concatenate three preprocessor tokens
#define CAT3_(a, b, c) a ## b ## c
#define CAT3(a, b, c) CAT3_(a, b, c)

//- Concatenate four preprocessor tokens
#define CAT4_(a, b, c, d) a ## b ## c ## d
#define CAT4(a, b, c, d) CAT4_(a, b, c, d)

//- Concatenate five preprocessor tokens
#define CAT5_(a, b, c, d, e) a ## b ## c ## d ## e
#define CAT5(a, b, c, d, e) CAT5_(a, b, c, d, e)

//- Generate an identifier unique within the file in which it is generated
#define FILE_UNIQUE(x) CAT(x, __LINE__)

//- Map 'name' to 'Name' via the predefined macro CAPITALIZE_name
#define CAPITALIZE(name) CAPITALIZE_##name

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
