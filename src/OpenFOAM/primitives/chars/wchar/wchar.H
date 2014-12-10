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

Primitive
    wchar_t

Description
    A wide-character and a pointer to a wide-character string.

SourceFiles
    wcharIO.C

SeeAlso
    http://en.wikipedia.org/wiki/UTF-8
    http://en.wikibooks.org/wiki/Unicode/Character_reference

\*---------------------------------------------------------------------------*/

#ifndef wchar_H
#define wchar_H

#include <cwchar>
#include <string>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Istream;
class Ostream;

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

//- Output wide character (Unicode) as UTF-8
Ostream& operator<<(Ostream&, const wchar_t);

//- Output wide character (Unicode) string as UTF-8
Ostream& operator<<(Ostream&, const wchar_t*);

//- Output wide character (Unicode) string as UTF-8
Ostream& operator<<(Ostream&, const std::wstring&);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
