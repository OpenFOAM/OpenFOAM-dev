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

Class
    Foam::Switch

Description
    A simple wrapper around bool so that it can be read as a word:
    true/false, on/off, yes/no, y/n, t/f, or none.

SourceFiles
    Switch.C
    SwitchIO.C

\*---------------------------------------------------------------------------*/

#ifndef Switch_H
#define Switch_H

#include "bool.H"
#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class Switch;
class dictionary;

Istream& operator>>(Istream&, Switch&);
Ostream& operator<<(Ostream&, const Switch&);


/*---------------------------------------------------------------------------*\
                           Class Switch Declaration
\*---------------------------------------------------------------------------*/

class Switch
{
    // Private data

        //- The logic and enumerated text representation stored as a single byte
        unsigned char switch_;

public:

    // Public data types

        // avoid issues with pre-processor defines
        #undef FALSE
        #undef TRUE
        #undef OFF
        #undef ON
        #undef NO
        #undef YES
        #undef NO_1
        #undef YES_1
        #undef FALSE_1
        #undef TRUE_1
        #undef NONE
        #undef PLACEHOLDER
        #undef INVALID

        //- The various text representations for a switch value.
        //  These also correspond to the entries in names.
        enum switchType
        {
            FALSE   = 0,  TRUE   = 1,
            OFF     = 2,  ON     = 3,
            NO      = 4,  YES    = 5,
            NO_1    = 6,  YES_1  = 7,
            FALSE_1 = 8,  TRUE_1 = 9,
            NONE    = 10, PLACEHOLDER = 11,
            INVALID
        };

    // Static data members

        //- The set of names corresponding to the switchType enumeration
        //  Includes an extra entry for "invalid".
        static const char* names[INVALID+1];


private:

    // Static Member Functions

        //- Return a switchType representation of a word
        static switchType asEnum(const std::string&, const bool allowInvalid);


public:

    // Constructors

        //- Construct null as false
        Switch()
        :
            switch_(Switch::FALSE)
        {}

        //- Construct from enumerated value
        Switch(const switchType sw)
        :
            switch_(sw)
        {}

        //- Construct from bool
        Switch(const bool b)
        :
            switch_(b ? Switch::TRUE : Switch::FALSE)
        {}

        //- Construct from integer values (treats integer as bool value)
        Switch(const int i)
        :
            switch_(i ? Switch::TRUE : Switch::FALSE)
        {}

        //- Construct from std::string, string, word
        //  Optionally allow bad words, and catch the error elsewhere
        Switch(const std::string& str, const bool allowInvalid=false)
        :
            switch_(asEnum(str, allowInvalid))
        {}

        //- Construct from character array
        //  Optionally allow bad words, and catch the error elsewhere
        Switch(const char* str, const bool allowInvalid=false)
        :
            switch_(asEnum(std::string(str), allowInvalid))
        {}

        //- Construct from Istream
        Switch(Istream& is);

        //- Construct from dictionary, supplying default value so that if the
        //  value is not found, it is added into the dictionary.
        static Switch lookupOrAddToDict
        (
            const word&,
            dictionary&,
            const Switch& defaultValue = false
        );


    // Member Functions

        //- Return true if the Switch has a valid value
        bool valid() const;

        //- Return a text representation of the Switch
        const char* asText() const;

        //- Update the value of the Switch if it is found in the dictionary
        bool readIfPresent(const word&, const dictionary&);


    // Member Operators

        //- Conversion to bool
        operator bool() const
        {
            return (switch_ & 0x1);
        }

        //- Assignment to enumerated value
        void operator=(const switchType sw)
        {
            switch_ = sw;
        }

        //- Assignment to bool
        void operator=(const bool b)
        {
            switch_ = (b ? Switch::TRUE : Switch::FALSE);
        }


    // IOstream Operators

        friend Istream& operator>>(Istream&, Switch&);
        friend Ostream& operator<<(Ostream&, const Switch&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
