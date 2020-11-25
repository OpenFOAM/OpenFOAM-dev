/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::wordRe

Description
    A wordRe is a word, but can also have a regular expression for matching
    words.

    By default the constructors will generally preserve the argument as a string
    literal and the assignment operators will use the wordRe::compOption::detect
    compOption to scan the string for regular expression meta characters and/or
    invalid word characters and react accordingly.

    The exceptions are when constructing/assigning from another
    Foam::wordRe (preserve the same type) or from a Foam::word (always
    literal).

Note
    If the string contents are changed - eg, by the operator+=() or by
    string::replace(), etc - it will be necessary to use compile() or
    recompile() to synchronize the regular expression.

SourceFiles
    wordRe.C

\*---------------------------------------------------------------------------*/

#ifndef wordRe_H
#define wordRe_H

#include "word.H"
#include "regExp.H"
#include "keyType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class wordRe;
class Istream;
class Ostream;

Istream& operator>>(Istream&, wordRe&);
Ostream& operator<<(Ostream&, const wordRe&);


/*---------------------------------------------------------------------------*\
                           Class wordRe Declaration
\*---------------------------------------------------------------------------*/

class wordRe
:
    public word
{
    // Private member data

        //- The regular expression
        mutable regExp re_;

public:

    // Static Data Members

        //- An empty wordRe
        static const wordRe null;


    // Public data types

        //- Enumeration with compile options
        //  Note that 'regexp' is implicit if 'noCase' is specified alone.
        enum class compOption
        {
            literal = 0, //!< treat as a string literal
            detect = 1,  //!< detect if the string contains meta-characters
            regExp = 2,  //!< treat as regular expression
            noCase  = 4, //!< ignore case in regular expression
            detectNoCase = detect | noCase,
            regExpNoCase = regExp | noCase
        };


        //- Is this a meta character?
        static inline bool meta(char);

        //- Test string for regular expression meta characters
        static inline bool isPattern(const string&);


    // Constructors

        //- Construct null
        inline wordRe();

        //- Copy constructor
        inline wordRe(const wordRe&);

        //- Construct from keyType
        inline wordRe(const keyType&);

        //- Construct from keyType
        inline wordRe(const keyType&, const compOption);

        //- Copy constructor of word
        inline wordRe(const word&);

        //- Copy constructor of character array
        //  Optionally specify how it should be treated.
        inline explicit wordRe
        (
            const char*,
            const compOption = compOption::literal
        );

        //- Copy constructor of string.
        //  Optionally specify how it should be treated.
        inline explicit wordRe
        (
            const string&,
            const compOption = compOption::literal
        );

        //- Copy constructor of std::string
        //  Optionally specify how it should be treated.
        inline explicit wordRe
        (
            const std::string&,
            const compOption = compOption::literal
        );

        //- Construct from Istream
        //  Words are treated as literals, strings with an auto-test
        wordRe(Istream&);


    // Member Functions

        // Access

            //- Should be treated as a match rather than a literal string?
            inline bool isPattern() const;


        // Infrastructure

            //- Compile the regular expression
            inline bool compile() const;

            //- Possibly compile the regular expression, with greater control
            inline bool compile(const compOption) const;

            //- Recompile an existing regular expression
            inline bool recompile() const;

            //- Frees precompiled regular expression, making wordRe a literal.
            //  Optionally strips invalid word characters
            inline void uncompile(const bool doStripInvalid = false) const;


        // Editing

            //- Copy string, auto-test for regular expression or other options
            inline void set
            (
                const std::string&,
                const compOption = compOption::detect
            );

            //- Copy string, auto-test for regular expression or other options
            inline void set
            (
                const char*,
                const compOption = compOption::detect
            );

            //- Clear string and precompiled regular expression
            inline void clear();


        // Searching

            //- Smart match as regular expression or as a string
            //  Optionally force a literal match only
            inline bool match
            (
                const std::string&,
                bool literalMatch = false
            ) const;


        // Miscellaneous

            //- Return a string with quoted meta-characters
            inline string quotemeta() const;

            //- Output some basic info
            Ostream& info(Ostream&) const;


    // Member Operators

        // Assignment

            //- Assign copy
            //  Always case sensitive
            inline void operator=(const wordRe&);

            //- Copy word, never a regular expression
            inline void operator=(const word&);

            //- Copy keyType, auto-test for regular expression
            //  Always case sensitive
            inline void operator=(const keyType&);

            //- Copy string, auto-test for regular expression
            //  Always case sensitive
            inline void operator=(const string&);

            //- Copy string, auto-test for regular expression
            //  Always case sensitive
            inline void operator=(const std::string&);

            //- Copy string, auto-test for regular expression
            //  Always case sensitive
            inline void operator=(const char*);


    // IOstream Operators

        friend Istream& operator>>(Istream&, wordRe&);
        friend Ostream& operator<<(Ostream&, const wordRe&);
};


inline int operator&
(
    const wordRe::compOption co1,
    const wordRe::compOption co2
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "wordReI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
