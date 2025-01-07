/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "LagrangianInjection.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class ... OtherFields>
void Foam::LagrangianInjection::filter
(
    barycentricField& coordinates,
    labelField& celli,
    labelField& facei,
    labelField& faceTrii,
    OtherFields& ... otherFields
)
{
    label iNew = 0;

    forAll(coordinates, iOld)
    {
        if (celli[iOld] != -1)
        {
            coordinates[iNew] = coordinates[iOld];
            celli[iNew] = celli[iOld];
            facei[iNew] = facei[iOld];
            faceTrii[iNew] = faceTrii[iOld];

            (void)std::initializer_list<nil>
            {(
                otherFields[iNew] = otherFields[iOld],
                nil()
            ) ... };

            ++ iNew;
        }
    }

    coordinates.resize(iNew);
    celli.resize(iNew);
    facei.resize(iNew);
    faceTrii.resize(iNew);

    (void)std::initializer_list<nil>
    {(
        otherFields.resize(iNew),
        nil()
    ) ... };
}


// ************************************************************************* //
