/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "surfacePatchIOList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(surfacePatchIOList, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from IOObject
Foam::surfacePatchIOList::surfacePatchIOList
(
    const IOobject& io
)
:
    surfacePatchList(),
    regIOobject(io)
{
    Foam::string functionName =
        "surfacePatchIOList::surfacePatchIOList"
        "(const IOobject& io)";


    if
    (
        readOpt() == IOobject::MUST_READ
     || readOpt() == IOobject::MUST_READ_IF_MODIFIED
    )
    {
        if (readOpt() == IOobject::MUST_READ_IF_MODIFIED)
        {
            WarningInFunction
                << "Specified IOobject::MUST_READ_IF_MODIFIED but class"
                << " does not support automatic rereading."
                << endl;
        }


        surfacePatchList& patches = *this;

        // read polyPatchList
        Istream& is = readStream(typeName);

        PtrList<entry> patchEntries(is);
        patches.setSize(patchEntries.size());

        label facei = 0;

        forAll(patches, patchi)
        {
            const dictionary& dict = patchEntries[patchi].dict();

            label patchSize = dict.lookup<label>("nFaces");
            label startFacei = dict.lookup<label>("startFace");

            patches[patchi] =
                surfacePatch
                (
                    word(dict.lookup("geometricType")),
                    patchEntries[patchi].keyword(),
                    patchSize,
                    startFacei,
                    patchi
                );


            if (startFacei != facei)
            {
                FatalErrorInFunction
                    << "Patches are not ordered. Start of patch " << patchi
                    << " does not correspond to sum of preceding patches."
                    << endl
                    << "while reading " << io.objectPath()
                    << exit(FatalError);
            }

            facei += patchSize;
        }

        // Check state of IOstream
        is.check(functionName.c_str());

        close();
    }
}

// Construct from IOObject
Foam::surfacePatchIOList::surfacePatchIOList
(
    const IOobject& io,
    const surfacePatchList& patches
)
:
    surfacePatchList(patches),
    regIOobject(io)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfacePatchIOList::~surfacePatchIOList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// writeData member function required by regIOobject
bool Foam::surfacePatchIOList::writeData(Ostream& os) const
{
    os  << *this;
    return os.good();
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const surfacePatchIOList& patches)
{
    os  << patches.size() << nl << token::BEGIN_LIST;

    forAll(patches, patchi)
    {
        patches[patchi].writeDict(os);
    }

    os  << token::END_LIST;

    return os;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// ************************************************************************* //
