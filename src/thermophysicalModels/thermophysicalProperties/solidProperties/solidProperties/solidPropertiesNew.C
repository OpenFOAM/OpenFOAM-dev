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

#include "solidProperties.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::solidProperties> Foam::solidProperties::New
(
    const word& name
)
{
    if (debug)
    {
        InfoInFunction << "Constructing solidProperties" << endl;
    }

    ConstructorTable::iterator cstrIter = ConstructorTablePtr_->find(name);

    if (cstrIter == ConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown solidProperties type "
            << name << nl << nl
            << "Valid solidProperties types are:" << nl
            << ConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<solidProperties>(cstrIter()());
}


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

    if (dict.found("defaultCoeffs"))
    {
        // Backward-compatibility

        if (Switch(dict.lookup("defaultCoeffs")))
        {
            return New(solidType);
        }
        else
        {
            return autoPtr<solidProperties>
            (
                new solidProperties(dict.optionalSubDict(solidType + "Coeffs"))
            );
        }
    }
    else
    {
        dictionaryConstructorTable::iterator cstrIter =
            dictionaryConstructorTablePtr_->find(solidType);

        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            FatalErrorInFunction
                << "Unknown solidProperties type "
                << solidType << nl << nl
                << "Valid solidProperties types are:" << nl
                << dictionaryConstructorTablePtr_->sortedToc()
                << exit(FatalError);
        }

        return autoPtr<solidProperties>(cstrIter()(dict));
    }
}


// ************************************************************************* //
