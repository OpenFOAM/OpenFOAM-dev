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

#include "ParticleForce.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::autoPtr<Foam::ParticleForce<CloudType>>
Foam::ParticleForce<CloudType>::New
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& name
)
{
    word forceType = name;
    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(forceType);

    if (cstrIter == dictionaryConstructorTablePtr_->end() && dict.found("type"))
    {
        forceType = dict.lookup<word>("type");
        cstrIter = dictionaryConstructorTablePtr_->find(forceType);
    }

    Info<< "    Selecting particle force " << forceType << endl;

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown particle force type "
            << forceType
            << ", constructor not in hash table" << nl << nl
            << "    Valid particle force types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<ParticleForce<CloudType>>
    (
        cstrIter()
        (
            owner,
            mesh,
            dict
        )
    );
}


// ************************************************************************* //
