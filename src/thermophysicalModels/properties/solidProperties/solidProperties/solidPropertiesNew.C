/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "solidProperties.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::solidProperties> Foam::solidProperties::New
(
    const dictionary& dict
)
{
    if (debug)
    {
        InfoInFunction << "Constructing solid" << endl;
    }

    const word solidType(dict.dictName());
    const Switch defaultCoeffs(dict.lookup("defaultCoeffs"));

    if (defaultCoeffs)
    {
        ConstructorTable::iterator cstrIter =
            ConstructorTablePtr_->find(solidType);

        if (cstrIter == ConstructorTablePtr_->end())
        {
            FatalErrorInFunction
                << "Unknown solidProperties type " << solidType << nl << nl
                << "Valid solidProperties types are :" << endl
                << ConstructorTablePtr_->sortedToc()
                << exit(FatalError);
        }

        return autoPtr<solidProperties>(cstrIter()());
    }
    else
    {
        return autoPtr<solidProperties>
        (
            new solidProperties
            (
                dict.subDict(solidType + "Coeffs")
            )
        );
    }
}


// ************************************************************************* //
