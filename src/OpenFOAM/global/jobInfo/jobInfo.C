/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#include "jobInfo.H"
#include "OSspecific.H"
#include "clock.H"
#include "OFstream.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

bool Foam::jobInfo::writeJobControl
(
    Foam::debug::infoSwitch("writeJobControl", 0)
);

bool Foam::jobInfo::writeJobInfo
(
    Foam::debug::infoSwitch("writeJobInfo", 0)
);

Foam::jobInfo Foam::jobInfo_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::jobInfo::jobInfo()
:
    runningJobPath_(),
    finishedJobPath_(),
    cpuTime_()
{
    name() = "jobInfo";

    if (Pstream::master())
    {
        if (writeJobControl)
        {
            string baseDir = getEnv("FOAM_JOB_DIR");
            string jobFile = hostName() + '.' + Foam::name(pid());

            fileName runningDir(baseDir/"runningJobs");
            fileName finishedDir(baseDir/"finishedJobs");

            runningJobPath_  = runningDir/jobFile;
            finishedJobPath_ = finishedDir/jobFile;

            if (baseDir.empty())
            {
                FatalErrorInFunction
                    << "Cannot get jobInfo directory $FOAM_JOB_DIR"
                    << Foam::exit(FatalError);
            }

            if (!isDir(runningDir) && !mkDir(runningDir))
            {
                FatalErrorInFunction
                    << "Cannot make jobInfo directory " << runningDir
                    << Foam::exit(FatalError);
            }

            if (!isDir(finishedDir) && !mkDir(finishedDir))
            {
                FatalErrorInFunction
                    << "Cannot make jobInfo directory " << finishedDir
                    << Foam::exit(FatalError);
            }

            writeJobInfo = true;
        }
    }

    constructed = true;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::jobInfo::~jobInfo()
{
    if (writeJobInfo && constructed && Pstream::master())
    {
        mv(runningJobPath_, finishedJobPath_);
    }

    constructed = false;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::jobInfo::write(Ostream& os) const
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


void Foam::jobInfo::write
(
    const word& executable,
    const fileName& casePath
) const
{
    if (writeJobInfo && Pstream::master())
    {
        if (!writeJobControl)
        {
            const fileName jobInfoPath(casePath/"jobInfo");

            if (!isDir(jobInfoPath) && !mkDir(jobInfoPath))
            {
                FatalErrorInFunction
                    << "Cannot make jobInfo directory " << jobInfoPath
                    << Foam::exit(FatalError);
            }

            const word jobFile = executable + '.' + Foam::name(pid());

            runningJobPath_ = jobInfoPath/jobFile;
            finishedJobPath_ = jobInfoPath/jobFile;
        }

        if (!write(OFstream(runningJobPath_)()))
        {
            FatalErrorInFunction
                << "Failed to write to jobInfo file "
                << runningJobPath_
                << Foam::exit(FatalError);
        }
    }
}


void Foam::jobInfo::end(const word& terminationType)
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


void Foam::jobInfo::end()
{
    end("normal");
}


void Foam::jobInfo::exit()
{
    end("exit");
}


void Foam::jobInfo::abort()
{
    end("abort");
}


void Foam::jobInfo::signalEnd() const
{
    if (writeJobControl && constructed && Pstream::master())
    {
        mv(runningJobPath_, finishedJobPath_);
    }

    constructed = false;
}


// ************************************************************************* //
