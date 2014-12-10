/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

Foam::autoPtr<Foam::solidProperties> Foam::solidProperties::New(Istream& is)
{
    if (debug)
    {
        Info<< "solidProperties::New(Istream&): constructing solid"
            << endl;
    }

    const word solidType(is);
    const word coeffs(is);

    if (coeffs == "defaultCoeffs")
    {
        ConstructorTable::iterator cstrIter =
            ConstructorTablePtr_->find(solidType);

        if (cstrIter == ConstructorTablePtr_->end())
        {
            FatalErrorIn("solidProperties::New(Istream&)")
                << "Unknown solidProperties type " << solidType << nl << nl
                << "Valid solidProperties types are :" << endl
                << ConstructorTablePtr_->sortedToc()
                << exit(FatalError);
        }

        return autoPtr<solidProperties>(cstrIter()());
    }
    else if (coeffs == "coeffs")
    {
        return autoPtr<solidProperties>(new solidProperties(is));
    }
    else
    {
        FatalErrorIn("solidProperties::New(Istream&)")
            << "solidProperties type " << solidType
            << ", option " << coeffs << " given"
            << ", should be coeffs or defaultCoeffs"
            << exit(FatalError);

        return autoPtr<solidProperties>(NULL);
    }
}


Foam::autoPtr<Foam::solidProperties> Foam::solidProperties::New
(
    const dictionary& dict
)
{
    if (debug)
    {
        Info<< "solidProperties::New(const dictionary&): constructing solid"
            << endl;
    }

    const word solidType(dict.dictName());
    const Switch defaultCoeffs(dict.lookup("defaultCoeffs"));

    if (defaultCoeffs)
    {
        ConstructorTable::iterator cstrIter =
            ConstructorTablePtr_->find(solidType);

        if (cstrIter == ConstructorTablePtr_->end())
        {
            FatalErrorIn("solidProperties::New(const dictionary&)")
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
