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

\*---------------------------------------------------------------------------*/

#include "JobInfo.H"
#include "OSspecific.H"
#include "clock.H"
#include "OFstream.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

bool Foam::JobInfo::writeJobInfo(Foam::debug::infoSwitch("writeJobInfo", 0));
Foam::JobInfo Foam::jobInfo;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
Foam::JobInfo::JobInfo()
:
    runningJobPath_(),
    finishedJobPath_(),
    cpuTime_()
{
    name() = "JobInfo";

    if (writeJobInfo && Pstream::master())
    {
        string baseDir = getEnv("FOAM_JOB_DIR");
        string jobFile = hostName() + '.' + Foam::name(pid());

        fileName runningDir(baseDir/"runningJobs");
        fileName finishedDir(baseDir/"finishedJobs");

        runningJobPath_  = runningDir/jobFile;
        finishedJobPath_ = finishedDir/jobFile;

        if (baseDir.empty())
        {
            FatalErrorIn("JobInfo::JobInfo()")
                << "Cannot get JobInfo directory $FOAM_JOB_DIR"
                << Foam::exit(FatalError);
        }

        if (!isDir(runningDir) && !mkDir(runningDir))
        {
            FatalErrorIn("JobInfo::JobInfo()")
                << "Cannot make JobInfo directory " << runningDir
                << Foam::exit(FatalError);
        }

        if (!isDir(finishedDir) && !mkDir(finishedDir))
        {
            FatalErrorIn("JobInfo::JobInfo()")
                << "Cannot make JobInfo directory " << finishedDir
                << Foam::exit(FatalError);
        }
    }

    constructed = true;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::JobInfo::~JobInfo()
{
    if (writeJobInfo && constructed && Pstream::master())
    {
        mv(runningJobPath_, finishedJobPath_);
    }

    constructed = false;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::JobInfo::write(Ostream& os) const
{
    if (writeJobInfo && Pstream::master())
    {
        if (os.good())
        {
            dictionary::write(os, false);
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return true;
    }
}


void Foam::JobInfo::write() const
{
    if (writeJobInfo && Pstream::master())
    {
        if (!write(OFstream(runningJobPath_)()))
        {
            FatalErrorIn("JobInfo::write() const")
                << "Failed to write to JobInfo file "
                << runningJobPath_
                << Foam::exit(FatalError);
        }
    }
}


void Foam::JobInfo::end(const word& terminationType)
{
    if (writeJobInfo && constructed && Pstream::master())
    {
        add("cpuTime", cpuTime_.elapsedCpuTime());
        add("endDate", clock::date());
        add("endTime", clock::clockTime());

        if (!found("termination"))
        {
            add("termination", terminationType);
        }

        rm(runningJobPath_);
        write(OFstream(finishedJobPath_)());
    }

    constructed = false;
}


void Foam::JobInfo::end()
{
    end("normal");
}


void Foam::JobInfo::exit()
{
    end("exit");
}


void Foam::JobInfo::abort()
{
    end("abort");
}


void Foam::JobInfo::signalEnd() const
{
    if (writeJobInfo && constructed && Pstream::master())
    {
        mv(runningJobPath_, finishedJobPath_);
    }

    constructed = false;
}


// ************************************************************************* //
