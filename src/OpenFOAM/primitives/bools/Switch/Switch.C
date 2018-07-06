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

\*---------------------------------------------------------------------------*/

#include "Switch.H"
#include "error.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* Foam::Switch::names[Foam::Switch::INVALID+1] =
{
    "false", "true",
    "off",   "on",
    "no",    "yes",
    "n",     "y",
    "f",     "t",
    "none",  "true",  // Is there a reasonable counterpart to "none"?
    "invalid"
};


// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

Foam::Switch::switchType Foam::Switch::asEnum
(
    const std::string& str,
    const bool allowInvalid
)
{
    for (int sw = 0; sw < Switch::INVALID; ++sw)
    {
        if (str == names[sw])
        {
            // handle aliases
            switch (sw)
            {
                case Switch::NO_1:
                case Switch::NONE:
                {
                    return Switch::NO;
                    break;
                }

                case Switch::YES_1:
                {
                    return Switch::YES;
                    break;
                }

                case Switch::FALSE_1:
                {
                    return Switch::FALSE;
                    break;
                }

                case Switch::TRUE_1:
                {
                    return Switch::TRUE;
                    break;
                }

                default:
                {
                    return switchType(sw);
                    break;
                }
            }
        }
    }

    if (!allowInvalid)
    {
        FatalErrorInFunction
            << "unknown switch word " << str << nl
            << abort(FatalError);
    }

    return Switch::INVALID;
}


Foam::Switch Foam::Switch::lookupOrAddToDict
(
    const word& name,
    dictionary& dict,
    const Switch& defaultValue
)
{
    return dict.lookupOrAddDefault<Switch>(name, defaultValue);
}


// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * //

bool Foam::Switch::valid() const
{
    return switch_ <= Switch::NONE;
}


const char* Foam::Switch::asText() const
{
    return names[switch_];
}


bool Foam::Switch::readIfPresent(const word& name, const dictionary& dict)
{
    return dict.readIfPresent<Switch>(name, *this);
}


// ************************************************************************* //
