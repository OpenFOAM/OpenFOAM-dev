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

#include "zoneGeneratorList.H"
#include "lookup.H"
#include "polyMesh.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::zoneGeneratorList::read
(
    const dictionary& dict,
    const bool topLevel
)
{
    setSize(dict.size());

    label i = 0;

    forAllConstIter(dictionary, dict, iter)
    {
        const word& name = iter().keyword();

        if (iter().isDict())
        {
            if (iter().dict().found("type"))
            {
                this->set
                (
                    i,
                    name,
                    zoneGenerator::New(name, mesh_, iter().dict())
                );
            }
            else
            {
                if (topLevel)
                {
                    // If the type is not specified
                    // assume the zoneGenerator type is the same as the name
                    this->set
                    (
                        i,
                        name,
                        zoneGenerator::New
                        (
                            name,
                            mesh_,
                            dictionary
                            (
                                iter().dict(),
                                dictionary::entries
                                (
                                    primitiveEntry
                                    (
                                        "type",
                                        name,
                                        iter().startLineNumber()
                                    )
                                )
                            )
                        )
                    );
                }
                else
                {
                    // If the type is not specified assume the it is a lookup
                    this->set
                    (
                        i,
                        name,
                        new zoneGenerators::lookup(name, mesh_, iter().dict())
                    );
                }
            }

            i++;
        }
        else if (!iter().stream().size())
        {
            if (topLevel)
            {
                // If an empty keyword is present assume it is a
                // zoneGenerator type and add to the list
                this->set
                (
                    i,
                    name,
                    zoneGenerator::New
                    (
                        name,
                        mesh_,
                        dictionary
                        (
                            name,
                            dict,
                            dictionary::entries
                            (
                                primitiveEntry
                                (
                                    "type",
                                    name,
                                    iter().startLineNumber()
                                )
                            )
                        )
                    )
                );
            }
            else
            {
                // If an empty keyword is present assume it is a zone name
                // and add a zone lookup
                this->set
                (
                    i,
                    name,
                    new zoneGenerators::lookup
                    (
                        name,
                        mesh_,
                        dictionary
                        (
                            name,
                            dict,
                            dictionary::entries
                            (
                                primitiveEntry
                                (
                                    "type",
                                    zoneGenerators::lookup::typeName,
                                    iter().startLineNumber()
                                )
                            )
                        )
                    )
                );
            }

            i++;
        }
    }

    setSize(i);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zoneGeneratorList::zoneGeneratorList(const polyMesh& mesh)
:
    PtrListDictionary<zoneGenerator>(0),
    mesh_(mesh)
{}


Foam::zoneGeneratorList::zoneGeneratorList
(
    const polyMesh& mesh,
    const dictionary& dict,
    const bool topLevel
)
:
    zoneGeneratorList(mesh)
{
    read(dict, topLevel);
}


Foam::zoneGeneratorList::zoneGeneratorList
(
    const polyMesh& mesh,
    const dictionary& dict,
    const word& zoneName,
    const word& zonesName
)
:
    zoneGeneratorList(mesh)
{
    if (dict.found(zoneName) && dict.found(zonesName))
    {
        FatalIOErrorInFunction(dict)
            << "Both " << zoneName << " and " << zonesName << " specified"
            << exit(FatalIOError);
    }

    if (dict.found(zoneName))
    {
        if (dict.isDict(zoneName))
        {
            append
            (
                zoneName,
                zoneGenerator::New(zoneName, mesh_, dict.subDict(zoneName))
            );
        }
        else
        {
            ITstream& zoneStream(dict.lookup(zoneName));
            const word name(zoneStream);
            append
            (
                name,
                new zoneGenerators::lookup
                (
                    name,
                    mesh,
                    dictionary
                    (
                        name,
                        dict,
                        dictionary::entries
                        (
                            primitiveEntry
                            (
                                "type",
                                zoneGenerators::lookup::typeName,
                                zoneStream.lineNumber()
                            )
                        )
                    )
                )
            );
        }
    }
    else if (dict.found(zonesName))
    {
        read(dict.subDict(zonesName));
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::zoneGeneratorList::moveUpdate() const
{
    const PtrListDictionary<zoneGenerator>& zgList(*this);

    forAll(zgList, i)
    {
        if (zgList[i].moveUpdate_)
        {
            return true;
        }
    }

    return false;
}


void Foam::zoneGeneratorList::generate()
{
    PtrListDictionary<zoneGenerator>& zgList(*this);

    forAll(zgList, i)
    {
        zgList[i].generate().store();
    }
}


void Foam::zoneGeneratorList::movePoints()
{
    PtrListDictionary<zoneGenerator>& zgList(*this);

    forAll(zgList, i)
    {
        zgList[i].movePoints().store();
    }
}


// ************************************************************************* //
