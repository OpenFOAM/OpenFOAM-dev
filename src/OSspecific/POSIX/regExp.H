/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::regExp

Description
    Wrapper around POSIX extended regular expressions.

    The PCRE '(?i)' extension is provided to compile the regular expression
    as being case-insensitive.

See also
    The manpage regex(7) for more information about POSIX regular expressions.
    These differ somewhat from \c Perl and \c sed regular expressions.

SourceFiles
    regExp.C

\*---------------------------------------------------------------------------*/

#ifndef regExp_H
#define regExp_H

#include <regex.h>
#include <string>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class string;
template<class T> class List;


/*---------------------------------------------------------------------------*\
                           Class regExp Declaration
\*---------------------------------------------------------------------------*/

class regExp
{
    // Private Data

        //- Precompiled regular expression
        mutable regex_t* preg_;


    // Private Member Functions

        //- Return true if it matches and sets the sub-groups matched.
        //  Templated to support both std::string and Foam::string
        template<class StringType>
        bool matchGrouping
        (
            const std::string&,
            List<StringType>& groups
        ) const;


public:

    // Static Member Functions

        //- Is character a regular expression meta-character?
        //  any character: '.' \n
        //  quantifiers: '*', '+', '?' \n
        //  grouping: '(', '|', ')' \n
        //  range: '[', ']' \n
        //
        //  Don't bother checking for '{digit}' bounds
        inline static bool meta(char c)
        {
            return
            (
                (c == '.')                           // any character
             || (c == '*' || c == '+' || c == '?')   // quantifiers
             || (c == '(' || c == ')' || c == '|')   // grouping/branching
             || (c == '[' || c == ']')               // range
            );
        }


    // Constructors

        //- Construct null
        regExp();

        //- Construct from character array, optionally ignoring case
        regExp(const char*, const bool ignoreCase=false);

        //- Construct from std::string (or string), optionally ignoring case
        regExp(const std::string&, const bool ignoreCase=false);

        //- Disallow default bitwise copy construction
        regExp(const regExp&) = delete;


    //- Destructor
    ~regExp();


    // Member Functions

        // Access

            //- Return true if a precompiled expression does not exist
            inline bool empty() const
            {
                return !preg_;
            }

            //- Does a precompiled expression exist?
            inline bool exists() const
            {
                return preg_ ? true : false;
            }

            //- Return the number of (groups)
            inline int ngroups() const
            {
                return int(preg_ ? preg_->re_nsub : 0);
            }


        // Editing

            //- Compile pattern into a regular expression,
            //  optionally ignoring case
            void set(const char*, const bool ignoreCase=false) const;

            //- Compile pattern into a regular expression,
            //  optionally ignoring case
            void set(const std::string&, const bool ignoreCase=false) const;

            //- Release precompiled expression.
            //  Returns true if precompiled expression existed before clear
            bool clear() const;


        // Searching

            //- Find position within string.
            //  Returns the index where it begins or string::npos if not found
            std::string::size_type find(const std::string& str) const;

            //- Return true if it matches the entire string
            //  The begin-of-line (^) and end-of-line ($) anchors are implicit
            bool match(const std::string&) const;

            //- Return true if it matches and sets the sub-groups matched
            //  The begin-of-line (^) and end-of-line ($) anchors are implicit
            bool match(const std::string&, List<std::string>& groups) const;

            //- Return true if it matches and sets the sub-groups matched
            //  The begin-of-line (^) and end-of-line ($) anchors are implicit
            bool match(const std::string&, List<string>& groups) const;

            //- Return true if the regex was found within string
            bool search(const std::string& str) const
            {
                return std::string::npos != find(str);
            }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const regExp&) = delete;

        //- Assign and compile pattern from a character array
        //  Always case sensitive
        void operator=(const char*);

        //- Assign and compile pattern from string
        //  Always case sensitive
        void operator=(const std::string&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
