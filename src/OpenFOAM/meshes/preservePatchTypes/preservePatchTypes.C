/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "preservePatchTypes.H"
#include "polyBoundaryMeshEntries.H"
#include "polyBoundaryMesh.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void Foam::preservePatchTypes
(
    const objectRegistry& obr,
    const word& meshInstance,
    const fileName& meshDir,
    const wordList& patchNames,
    PtrList<dictionary>& patchDicts,
    const word& defaultFacesName,
    word& defaultFacesType
)
{
    patchDicts.setSize(patchNames.size());

    dictionary patchDictionary;

    // Read boundary file as single dictionary
    {
        typeIOobject<polyBoundaryMeshEntries> patchEntriesHeader
        (
            "boundary",
            meshInstance,
            meshDir,
            obr,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );

        if (patchEntriesHeader.headerOk())
        {
            // Create a list of entries from the boundary file.
            polyBoundaryMeshEntries patchEntries(patchEntriesHeader);

            forAll(patchEntries, patchi)
            {
                patchDictionary.add(patchEntries[patchi]);
            }
        }
    }

    forAll(patchNames, patchi)
    {
        if (patchDictionary.found(patchNames[patchi]))
        {
            const dictionary& patchDict =
                patchDictionary.subDict(patchNames[patchi]);

            patchDicts.set(patchi, patchDict.clone());
            patchDicts[patchi].remove("nFaces");
            patchDicts[patchi].remove("startFace");
        }
    }

    if (patchDictionary.found(defaultFacesName))
    {
        const dictionary& patchDict =
            patchDictionary.subDict(defaultFacesName);

        patchDict.readIfPresent("geometricType", defaultFacesType);
    }

    Info<< nl << "Default patch type set to " << defaultFacesType << endl;
}


// ************************************************************************* //
