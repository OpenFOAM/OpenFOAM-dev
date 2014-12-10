/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::word Foam::Time::findInstance
(
    const fileName& dir,
    const word& name,
    const IOobject::readOption rOpt,
    const word& stopInstance
) const
{
    // Note: if name is empty, just check the directory itself


    const fileName tPath(path());
    const fileName dirPath(tPath/timeName()/dir);

    // check the current time directory
    if
    (
        name.empty()
      ? isDir(dirPath)
      :
        (
            isFile(dirPath/name)
         && IOobject(name, timeName(), dir, *this).headerOk()
        )
    )
    {
        if (debug)
        {
            Info<< "Time::findInstance"
                "(const fileName&, const word&"
                ", const IOobject::readOption, const word&)"
                << " : found \"" << name
                << "\" in " << timeName()/dir
                << endl;
        }

        return timeName();
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
        if
        (
            name.empty()
          ? isDir(tPath/ts[instanceI].name()/dir)
          :
            (
                isFile(tPath/ts[instanceI].name()/dir/name)
             && IOobject(name, ts[instanceI].name(), dir, *this).headerOk()
            )
        )
        {
            if (debug)
            {
                Info<< "Time::findInstance"
                    "(const fileName&, const word&"
                    ", const IOobject::readOption, const word&)"
                    << " : found \"" << name
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
                Info<< "Time::findInstance"
                    "(const fileName&, const word&"
                    ", const IOobject::readOption, const word&)"
                    << " : hit stopInstance " << stopInstance
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
                    FatalErrorIn
                    (
                        "Time::findInstance"
                        "(const fileName&, const word&"
                        ", const IOobject::readOption, const word&)"
                    )   << "Cannot find directory "
                        << dir << " in times " << timeName()
                        << " down to " << stopInstance
                        << exit(FatalError);
                }
                else
                {
                    FatalErrorIn
                    (
                        "Time::findInstance"
                        "(const fileName&, const word&"
                        ", const IOobject::readOption, const word&)"
                    )   << "Cannot find file \"" << name << "\" in directory "
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

    if
    (
        name.empty()
      ? isDir(tPath/constant()/dir)
      :
        (
            isFile(tPath/constant()/dir/name)
         && IOobject(name, constant(), dir, *this).headerOk()
        )
    )
    {
        if (debug)
        {
            Info<< "Time::findInstance"
                "(const fileName&, const word&"
                ", const IOobject::readOption, const word&)"
                << " : found \"" << name
                << "\" in " << constant()/dir
                << endl;
        }

        return constant();
    }

    if (rOpt == IOobject::MUST_READ || rOpt == IOobject::MUST_READ_IF_MODIFIED)
    {
        FatalErrorIn
        (
            "Time::findInstance"
            "(const fileName&, const word&"
            ", const IOobject::readOption, const word&)"
        )   << "Cannot find file \"" << name << "\" in directory "
            << dir << " in times " << timeName()
            << " down to " << constant()
            << exit(FatalError);
    }

    return constant();
}


// ************************************************************************* //
