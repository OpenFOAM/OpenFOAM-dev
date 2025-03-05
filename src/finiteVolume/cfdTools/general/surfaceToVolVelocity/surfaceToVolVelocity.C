/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2025 OpenFOAM Foundation
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

#include "surfaceToVolVelocity.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

const Foam::volVectorField& Foam::surfaceToVolVelocity
(
    const surfaceVectorField& Uf
)
{
    if
    (
        Uf.dimensions() != dimVelocity
     && Uf.dimensions() != dimDensity*dimVelocity
     && Uf.dimensions() != dimVelocity/dimTime
     && Uf.dimensions() != dimDensity*dimVelocity/dimTime
    )
    {
        return volVectorField::null();
    }

    const word& surfaceName = Uf.name();

    // Split the name based on non-alphanumeric (and then some) separators
    DynamicList<word> parts;
    DynamicList<char> separators;
    {
        string::size_type i0 = 0;

        for (string::size_type i = 0; i < surfaceName.size(); ++ i)
        {
            const char c = surfaceName[i];

            if (isalnum(c) || c == '_' || c == '.') continue;

            parts.append(surfaceName(i0, i - i0));
            separators.append(c);

            i0 = i + 1;
        }

        parts.append(surfaceName(i0, surfaceName.size() - i0));
    }

    // Convert every part to its equivalent vol name
    forAll(parts, parti)
    {
        word& part = parts[parti];

        // Find where the groups start
        word::size_type i = part.find_first_of('.');
        if (i == word::npos) i = part.size();

        // Find where the old-time suffixes start
        while (i > 0 && part(i - 2, 2) == "_0") i -= 2;

        // We are back to the end of the core field name. If this ends with an
        // 'f' then this should be removed to form the equivalent vol name.
        if (i > 0 && part[i - 1] == 'f')
        {
            part = word(part(i - 1) + part(i, word::npos));
        }
    }

    // Put the parts back together and return
    word volName;
    forAll(separators, parti)
    {
        volName += parts[parti] + separators[parti];
    }
    volName += parts.last();

    // If no 'f' suffixes were removed then there is no corresponding vol name
    if (volName == surfaceName)
    {
        return volVectorField::null();
    }

    // Return the vol field if it can be found
    return
        Uf.mesh().foundObject<volVectorField>(volName)
      ? Uf.mesh().lookupObject<volVectorField>(volName)
      : volVectorField::null();
}


// ************************************************************************* //
