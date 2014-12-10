/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Global
    writeMeshObject

Description
    Write a mesh object in the format specified in controlDict.

\*---------------------------------------------------------------------------*/

#ifndef writeMeshObject_H
#define writeMeshObject_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class T>
inline bool writeMeshObject
(
    const word& name,
    const fileName& meshDir,
    Time& runTime,
    bool disableHeaderChecking = false
)
{
    IOobject io
    (
        name,
        runTime.timeName(),
        meshDir,
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    bool writeOk = false;

    if (io.headerOk())
    {
        Info<< "        Reading " << io.headerClassName()
            << " : " << name << endl;

        // Switch off type checking (for reading e.g. faceZones as
        // generic list of dictionaries).
        word oldTypeName;
        if (disableHeaderChecking)
        {
            oldTypeName = T::typeName;
            const_cast<word&>(T::typeName) = word::null;
        }

        T meshObject(io);

        if (disableHeaderChecking)
        {
            const_cast<word&>(T::typeName) = oldTypeName;
            // Fake type back to what was in field
            const_cast<word&>(meshObject.type()) = io.headerClassName();
        }

        Info<< "        Writing " << name << endl;
        writeOk = meshObject.regIOobject::write();
    }

    return writeOk;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
