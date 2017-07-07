/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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
    If "name" is empty: return the location of "directory"
    If "name" is not empty: return the location of "directory" containing the
    file "name".
    Used in reading mesh data.

\*---------------------------------------------------------------------------*/

#include "Time.H"
#include "IOobject.H"
#include "IOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool Foam::Time::exists(IOobject& io)
{
    // Generate filename for object
    fileName objPath(fileHandler().objectPath(io, word::null));

    // Test for either directory or a (valid) file & IOobject
    bool ok;
    if (io.name().empty())
    {
        ok = fileHandler().isDir(objPath);
    }
    else
    {
        ok =
            fileHandler().isFile(objPath)
         && io.typeHeaderOk<IOList<label>>(false);// object with local scope
    }

    if (!ok)
    {
        // Re-test with raw objectPath. This is for backwards
        // compatibility
        fileName originalPath(io.objectPath());
        if (originalPath != objPath)
        {
            // Test for either directory or a (valid) file & IOobject
            if (io.name().empty())
            {
                ok = fileHandler().isDir(originalPath);
            }
            else
            {
                ok =
                    fileHandler().isFile(originalPath)
                 && io.typeHeaderOk<IOList<label>>(false);
            }
        }
    }

    return ok;
}


Foam::word Foam::Time::findInstance
(
    const fileName& dir,
    const word& name,
    const IOobject::readOption rOpt,
    const word& stopInstance
) const
{
    // Note: - if name is empty, just check the directory itself
    //       - check both for isFile and headerOk since the latter does a
    //         filePath so searches for the file.
    //       - check for an object with local file scope (so no looking up in
    //         parent directory in case of parallel)

    {
        IOobject io
        (
            name,           // name might be empty!
            timeName(),
            dir,
            *this
        );

        if (exists(io))
        {
            if (debug)
            {
                InfoInFunction
                    << "Found exact match for \"" << name
                    << "\" in " << timeName()/dir
                    << endl;
            }

            return timeName();
        }
    }

    // Search back through the time directories to find the time
    // closest to and lower than current time

    instantList ts = times();
    label instanceI;

    for (instanceI = ts.size()-1; instanceI >= 0; --instanceI)
    {
        if (ts[instanceI].value() <= timeOutputValue())
        {
            break;
        }
    }

    // continue searching from here
    for (; instanceI >= 0; --instanceI)
    {
        IOobject io
        (
            name,           // name might be empty!
            ts[instanceI].name(),
            dir,
            *this
        );

        if (exists(io))
        {
            if (debug)
            {
                InfoInFunction
                    << "Found instance match for \"" << name
                    << "\" in " << ts[instanceI].name()/dir
                    << endl;
            }

            return ts[instanceI].name();
        }

        // Check if hit minimum instance
        if (ts[instanceI].name() == stopInstance)
        {
            if (debug)
            {
                //InfoInFunction
                Pout<< "findInstance : "
                    << "Hit stopInstance " << stopInstance
                    << endl;
            }

            if
            (
                rOpt == IOobject::MUST_READ
             || rOpt == IOobject::MUST_READ_IF_MODIFIED
            )
            {
                if (name.empty())
                {
                    FatalErrorInFunction
                        << "Cannot find directory "
                        << dir << " in times " << timeName()
                        << " down to " << stopInstance
                        << exit(FatalError);
                }
                else
                {
                    FatalErrorInFunction
                        << "Cannot find file \"" << name << "\" in directory "
                        << dir << " in times " << timeName()
                        << " down to " << stopInstance
                        << exit(FatalError);
                }
            }

            return ts[instanceI].name();
        }
    }


    // not in any of the time directories, try constant

    // Note. This needs to be a hard-coded constant, rather than the
    // constant function of the time, because the latter points to
    // the case constant directory in parallel cases

    IOobject io
    (
        name,
        constant(),
        dir,
        *this
    );

    if (exists(io))
    {
        if (debug)
        {
            InfoInFunction
                << "Found constant match for \"" << name
                << "\" in " << constant()/dir
                << endl;
        }

        return constant();
    }

    if (rOpt == IOobject::MUST_READ || rOpt == IOobject::MUST_READ_IF_MODIFIED)
    {
        FatalErrorInFunction
            << "Cannot find file \"" << name << "\" in directory "
            << dir << " in times " << timeName()
            << " down to " << constant()
            << exit(FatalError);
    }

    return constant();
}


// ************************************************************************* //
