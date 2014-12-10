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

\*---------------------------------------------------------------------------*/

#include "triSurface.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool triSurface::readSTL(const fileName& STLfileName)
{
    token firstToken;
    {
        IFstream str(STLfileName);

        if (!str.good())
        {
            return false;
        }
        firstToken = token(str);
    }

    if (firstToken.isWord())
    {
        word firstWord(firstToken.wordToken());

        for (size_t i = 0; i < firstWord.size(); i++)
        {
            firstWord[i] = std::toupper(firstWord[i]);
        }

        if (firstWord == "SOLID" || firstWord(5) == "SOLID")
        {
            return readSTLASCII(STLfileName);
        }
    }

    return readSTLBINARY(STLfileName);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
