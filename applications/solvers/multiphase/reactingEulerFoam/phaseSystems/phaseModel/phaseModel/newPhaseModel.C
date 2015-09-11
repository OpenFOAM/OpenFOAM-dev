/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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

#include "phaseModel.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::phaseModel> Foam::phaseModel::New
(
    const phaseSystem& fluid,
    const word& phaseName,
    const label index
)
{
    word phaseModelType(fluid.subDict(phaseName).lookup("type"));

    Info<< "Selecting phaseModel for "
        << phaseName << ": " << phaseModelType << endl;

    phaseSystemConstructorTable::iterator cstrIter =
        phaseSystemConstructorTablePtr_->find(phaseModelType);

    if (cstrIter == phaseSystemConstructorTablePtr_->end())
    {
        FatalErrorIn("phaseModel::New")
            << "Unknown phaseModelType type "
            << phaseModelType << endl << endl
            << "Valid phaseModel types are : " << endl
            << phaseSystemConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(fluid, phaseName, index);
}

// ************************************************************************* //
