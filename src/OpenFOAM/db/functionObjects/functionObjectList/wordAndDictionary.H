/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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
    Foam::wordAndDictionary

Description
    Tuple of a word and dictionary, used to read in per-field options for
    function objects in the following syntax:

    fields
    (
        p
        {
            option1     true;
            option2     false;
        }

        U

        T
        {
            option1     false;
        }
    );

    IO is like the tuple, except that there are no enclosing parentheses, and
    if the dictionary is empty then '{}' is omitted. The latter means that in
    the absence of any options the syntax becomes that of a wordList, which
    means it can be used for argument parsing for input of that form.

SourceFiles
    wordAndDictionary.C

\*---------------------------------------------------------------------------*/

#ifndef wordAndDictionary_H
#define wordAndDictionary_H

#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class wordAndDictionary;
Istream& operator>>(Istream&, wordAndDictionary&);
Ostream& operator<<(Ostream&, const wordAndDictionary&);


/*---------------------------------------------------------------------------*\
                       Class wordAndDictionary Declaration
\*---------------------------------------------------------------------------*/

class wordAndDictionary
:
    public Tuple2<word, dictionary>
{
public:

    // Constructors

        //- Inherit tuple constructors
        using Tuple2<word, dictionary>::Tuple2;

        //- Construct null
        wordAndDictionary();

        //- Construct from Istream
        wordAndDictionary(Istream&);


    //- Destructor
    ~wordAndDictionary();


    // IOstream Operators

        friend Istream& operator>>(Istream&, wordAndDictionary&);
        friend Ostream& operator<<(Ostream&, const wordAndDictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
