/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "IOobject.H"
#include "fileOperation.H"
#include "Istream.H"
#include "IOstreams.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::IOobject::typeHeaderOk(const bool checkType)
{
    bool ok = true;

    // Everyone check or just master
    bool masterOnly =
        typeGlobal<Type>()
     && (
            IOobject::fileModificationChecking == timeStampMaster
         || IOobject::fileModificationChecking == inotifyMaster
        );

    const fileOperation& fp = Foam::fileHandler();

    // Determine local status
    if (!masterOnly || Pstream::master())
    {
        fileName fName(typeFilePath<Type>(*this));

        ok = fp.readHeader(*this, fName, Type::typeName);
        if (ok && checkType && headerClassName_ != Type::typeName)
        {
            WarningInFunction
                << "unexpected class name " << headerClassName_
                << " expected " << Type::typeName
                << " when reading " << fName << endl;

            ok = false;
        }
    }

    // If masterOnly make sure all processors know about it
    if (masterOnly)
    {
        Pstream::scatter(ok);
    }

    return ok;
}


template<class Type>
void Foam::IOobject::warnNoRereading() const
{
    if (readOpt() == IOobject::MUST_READ_IF_MODIFIED)
    {
        WarningInFunction
            << Type::typeName << ' ' << name()
            << " constructed with IOobject::MUST_READ_IF_MODIFIED"
            " but " << Type::typeName
            << " does not support automatic rereading."
            << endl;
    }
}


// ************************************************************************* //
