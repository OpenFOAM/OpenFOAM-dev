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

#include "surfZoneIOList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(surfZoneIOList, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfZoneIOList::surfZoneIOList
(
    const IOobject& io
)
:
    surfZoneList(),
    regIOobject(io)
{
    Foam::string functionName =
        "surfZoneIOList::surfZoneIOList"
        "(const IOobject& io)";


    if
    (
        readOpt() == IOobject::MUST_READ
     || readOpt() == IOobject::MUST_READ_IF_MODIFIED
    )
    {
        surfZoneList& zones = *this;

        Istream& is = readStream(typeName);

        PtrList<entry> dictEntries(is);
        zones.setSize(dictEntries.size());

        label facei = 0;
        forAll(zones, zoneI)
        {
            const dictionary& dict = dictEntries[zoneI].dict();

            label zoneSize = readLabel(dict.lookup("nFaces"));
            label startFacei = readLabel(dict.lookup("startFace"));

            zones[zoneI] = surfZone
            (
                dictEntries[zoneI].keyword(),
                zoneSize,
                startFacei,
                zoneI
            );

            word geoType;
            if (dict.readIfPresent("geometricType", geoType))
            {
                zones[zoneI].geometricType() = geoType;
            }

            if (startFacei != facei)
            {
                FatalErrorInFunction
                    << "surfZones are not ordered. Start of zone " << zoneI
                    << " does not correspond to sum of preceding zones." << nl
                    << "while reading " << io.objectPath() << endl
                    << exit(FatalError);
            }

            facei += zoneSize;
        }

        // Check state of IOstream
        is.check(functionName.c_str());

        close();
    }
}


Foam::surfZoneIOList::surfZoneIOList
(
    const IOobject& io,
    const surfZoneList& zones
)
:
    surfZoneList(zones),
    regIOobject(io)
{}


Foam::surfZoneIOList::surfZoneIOList
(
    const IOobject& io,
    surfZoneList&& zones
)
:
    surfZoneList(move(zones)),
    regIOobject(io)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfZoneIOList::~surfZoneIOList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// writeData member function required by regIOobject
bool Foam::surfZoneIOList::writeData(Ostream& os) const
{
    os  << *this;
    return os.good();
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const surfZoneIOList& L)
{
    os  << L.size() << nl << token::BEGIN_LIST << incrIndent << nl;

    forAll(L, i)
    {
        L[i].writeDict(os);
    }

    os  << decrIndent << token::END_LIST;

    return os;
}


// ************************************************************************* //
