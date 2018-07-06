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

#include "polyPatch.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::polyPatch> Foam::polyPatch::New
(
    const word& patchType,
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm
)
{
    if (debug)
    {
        InfoInFunction << "Constructing polyPatch" << endl;
    }

    wordConstructorTable::iterator cstrIter =
        wordConstructorTablePtr_->find(patchType);

    if (cstrIter == wordConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown polyPatch type "
            << patchType << " for patch " << name << nl << nl
            << "Valid polyPatch types are :" << endl
            << wordConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<polyPatch>
    (
        cstrIter()
        (
            name,
            size,
            start,
            index,
            bm,
            patchType
        )
    );
}


Foam::autoPtr<Foam::polyPatch> Foam::polyPatch::New
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm
)
{
    if (debug)
    {
        InfoInFunction << "Constructing polyPatch" << endl;
    }

    word patchType(dict.lookup("type"));
    dict.readIfPresent("geometricType", patchType);

    return polyPatch::New(patchType, name, dict, index, bm);
}


Foam::autoPtr<Foam::polyPatch> Foam::polyPatch::New
(
    const word& patchType,
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm
)
{
    if (debug)
    {
        InfoInFunction << "Constructing polyPatch" << endl;
    }

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(patchType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        if (!disallowGenericPolyPatch)
        {
            cstrIter = dictionaryConstructorTablePtr_->find("genericPatch");
        }

        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            FatalIOErrorInFunction
            (
                dict
            )   << "Unknown polyPatch type "
                << patchType << " for patch " << name << nl << nl
                << "Valid polyPatch types are :" << endl
                << dictionaryConstructorTablePtr_->sortedToc()
                << exit(FatalIOError);
        }
    }

    return autoPtr<polyPatch>(cstrIter()(name, dict, index, bm, patchType));
}


// ************************************************************************* //
