/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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

Description
    Writes the header description of the File to the stream
    associated with the File.

\*---------------------------------------------------------------------------*/

#include "IOobject.H"
#include "objectRegistry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool Foam::IOobject::writeHeader
(
    Ostream& os,
    const IOstream::versionNumber version,
    const IOstream::streamFormat format,
    const word& type,
    const string& note,
    const fileName& location,
    const word& name
)
{
    if (!os.good())
    {
        InfoInFunction
            << "No stream open for write" << nl
            << os.info() << endl;

        return false;
    }

    writeBanner(os) << foamFile << "\n{\n";

    if (version != IOstream::currentVersion)
    {
        os  << "    version     " << version << ";\n";
    }

    os  << "    format      " << format << ";\n"
        << "    class       " << type << ";\n";

    if (note.size())
    {
        os  << "    note        " << note << ";\n";
    }

    if (location.size())
    {
        os  << "    location    " << location << ";\n";
    }

    os  << "    object      " << name << ";\n"
        << "}" << nl;

    writeDivider(os) << nl;

    return true;
}


bool Foam::IOobject::writeHeader
(
    Ostream& os,
    const dictionary& foamFileDict
)
{
    return writeHeader
    (
        os,
        foamFileDict.lookupOrDefault<IOstream::versionNumber>
        (
            "version",
            IOstream::currentVersion
        ),
        foamFileDict.found("format")
      ? IOstream::formatEnum(foamFileDict.lookup<word>("format"))
      : IOstream::ASCII,
        foamFileDict.lookupOrDefault<word>("class", dictionary::typeName),
        foamFileDict.lookupOrDefault<string>("note", string::null),
        foamFileDict.lookupOrDefault<fileName>("location", fileName::null),
        foamFileDict.lookupOrDefault<word>("name", os.name().name())
    );
}


bool Foam::IOobject::writeHeader(Ostream& os, const word& type) const
{
    return writeHeader
    (
        os,
        os.version(),
        os.format(),
        type,
        note(),
        instance()/db().dbDir()/local(),
        name()
    );
}


bool Foam::IOobject::writeHeader(Ostream& os) const
{
    return writeHeader(os, type());
}


// ************************************************************************* //
