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

#include "tetherPotentialList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::tetherPotentialList::readTetherPotentialDict
(
    const List<word>& siteIdList,
    const dictionary& tetherPotentialDict,
    const List<word>& tetherSiteIdList
)
{

    Info<< nl << "Building tether potentials." << endl;

    idMap_ = List<label>(siteIdList.size(), -1);

    label tetherMapIndex = 0;

    forAll(tetherSiteIdList, t)
    {
        word tetherPotentialName = tetherSiteIdList[t];

        label tetherId = findIndex(siteIdList, tetherPotentialName);

        if (tetherId == -1)
        {
            FatalErrorInFunction
                << nl
                << "No matching entry found in siteIdList for tether name "
                << tetherPotentialName
                << abort(FatalError);
        }
        else if (!tetherPotentialDict.found(tetherPotentialName))
        {
            FatalErrorInFunction
                << nl << "tether potential specification subDict "
                << tetherPotentialName << " not found"
                << abort(FatalError);
        }
        else
        {
            this->set
            (
                tetherMapIndex,
                tetherPotential::New
                (
                    tetherPotentialName,
                    tetherPotentialDict.subDict(tetherPotentialName)
                )
            );
        }

        idMap_[tetherId] = tetherMapIndex;

        tetherMapIndex++;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::tetherPotentialList::tetherPotentialList()
:
    PtrList<tetherPotential>(),
    idMap_()
{}


Foam::tetherPotentialList::tetherPotentialList
(
    const List<word>& siteIdList,
    const dictionary& tetherPotentialDict,
    const List<word>& tetherSiteIdList
)
:
    PtrList<tetherPotential>(),
    idMap_()
{
    buildPotentials(siteIdList, tetherPotentialDict, tetherSiteIdList);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::tetherPotentialList::~tetherPotentialList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::tetherPotentialList::buildPotentials
(
    const List<word>& siteIdList,
    const dictionary& tetherPotentialDict,
    const List<word>& tetherSiteIdList
)
{
    setSize(tetherSiteIdList.size());

    readTetherPotentialDict(siteIdList, tetherPotentialDict, tetherSiteIdList);
}


const Foam::tetherPotential& Foam::tetherPotentialList::tetherPotentialFunction
(
    const label a
) const
{
    return (*this)[tetherPotentialIndex(a)];
}


Foam::vector Foam::tetherPotentialList::force
(
    const label a,
    const vector rIT
) const
{
    return (*this)[tetherPotentialIndex(a)].force(rIT);
}


Foam::scalar Foam::tetherPotentialList::energy
(
    const label a,
    const vector rIT
) const
{
    return (*this)[tetherPotentialIndex(a)].energy(rIT);
}


// ************************************************************************* //
