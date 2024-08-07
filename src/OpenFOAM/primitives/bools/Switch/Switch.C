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

\*---------------------------------------------------------------------------*/

#include "Switch.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* Foam::Switch::names[nSwitchType] =
{
    "false",
    "true",
    "off",
    "on",
    "no",
    "yes",
    "n",
    "y",
    "f",
    "t",
    "none",
    "any",
    "invalid"
};


// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

Foam::Switch::switchType Foam::Switch::asEnum
(
    const std::string& str,
    const bool allowInvalid
)
{
    for (switchType sw=switchType::False; sw<switchType::invalid; ++sw)
    {
        if (str == names[toInt(sw)])
        {
            // Handle aliases
            switch (sw)
            {
                case switchType::n:
                case switchType::none:
                {
                    return switchType::no;
                    break;
                }

                case switchType::y:
                case switchType::any:
                {
                    return switchType::yes;
                    break;
                }

                case switchType::f:
                {
                    return switchType::False;
                    break;
                }

                case switchType::t:
                {
                    return switchType::True;
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

    return switchType::invalid;
}


// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * //

bool Foam::Switch::valid() const
{
    return switch_ <= switchType::none;
}


const char* Foam::Switch::asText() const
{
    return names[toInt(switch_)];
}


bool Foam::Switch::readIfPresent(const word& name, const dictionary& dict)
{
    return dict.readIfPresent<Switch>(name, *this);
}


// ************************************************************************* //
