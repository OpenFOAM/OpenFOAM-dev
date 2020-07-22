/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "error.H"
#include "integrationScheme.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::integrationScheme> Foam::integrationScheme::New
(
    const word& phiName,
    const dictionary& dict
)
{
    const word schemeName(dict.lookup(phiName));

    Info<< "Selecting " << phiName << " integration scheme "
        << schemeName << endl;

    typename wordConstructorTable::iterator cstrIter =
        wordConstructorTablePtr_->find(schemeName);

    if (cstrIter == wordConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown integration scheme type "
            << schemeName << nl << nl
            << "Valid integration scheme types are:" << nl
            << wordConstructorTablePtr_->sortedToc() << nl
            << exit(FatalError);
    }

    return autoPtr<integrationScheme>(cstrIter()());
}

// ************************************************************************* //
