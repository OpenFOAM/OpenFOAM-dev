/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "progressFunctionObject.H"
#include "Time.H"
#include "OFstream.H"
#include "addToRunTimeSelectionTable.H"
#include <iomanip>
#include <sstream>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(progress, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        progress,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::string Foam::functionObjects::progress::formatTime(const scalar seconds) const
{
    if (seconds < 0 || seconds > 86400 * 365)
    {
        return "--:--";
    }

    const label totalSec = label(seconds + 0.5);
    const label hours = totalSec / 3600;
    const label minutes = (totalSec % 3600) / 60;
    const label secs = totalSec % 60;

    std::ostringstream oss;
    if (hours > 0)
    {
        oss << hours << "h "
            << std::setw(2) << std::setfill('0') << minutes << "m "
            << std::setw(2) << std::setfill('0') << secs << "s";
    }
    else if (minutes > 0)
    {
        oss << minutes << "m "
            << std::setw(2) << std::setfill('0') << secs << "s";
    }
    else
    {
        oss << secs << "s";
    }

    return oss.str();
}


void Foam::functionObjects::progress::writeFileHeader(const label i)
{
    if (Pstream::master())
    {
        writeHeader(file(), "Simulation Progress");
        writeCommented(file(), "Time");

        file() << tab << "Percent"
               << tab << "ElapsedClockTime[s]"
               << tab << "SolvingRate[sim_s/wall_s]"
               << tab << "ETA[s]";

        file() << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::progress::progress
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObject(name, runTime, dict),
    logFiles(runTime, name),
    startTime_(runTime.startTime().value()),
    endTime_(runTime.endTime().value()),
    log_(true),
    writeJson_(true),
    barWidth_(25)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::progress::~progress()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::progress::read(const dictionary& dict)
{
    functionObject::read(dict);

    resetName(typeName);

    dict.readIfPresent("log", log_);
    dict.readIfPresent("writeJson", writeJson_);
    dict.readIfPresent("barWidth", barWidth_);

    startTime_ = time_.startTime().value();
    endTime_ = time_.endTime().value();

    return true;
}


bool Foam::functionObjects::progress::execute()
{
    return true;
}


bool Foam::functionObjects::progress::write()
{
    logFiles::write();

    const scalar curTime = time_.value();
    const scalar totalSim = max(VSMALL, endTime_ - startTime_);
    const scalar completedSim = max(0.0, curTime - startTime_);
    const scalar percent = min(100.0, max(0.0, (completedSim / totalSim) * 100.0));

    const scalar elapsedClock = scalar(time_.elapsedClockTime());
    const scalar rate = completedSim / max(1.0, elapsedClock);

    scalar etaSeconds = 0.0;
    string etaString = "Completed";

    if (curTime < endTime_ && rate > VSMALL)
    {
        etaSeconds = (endTime_ - curTime) / rate;
        etaString = formatTime(etaSeconds);
    }

    if (Pstream::master())
    {
        // 1. Write tabular log file data
        writeTime(file());
        file() << tab << percent
               << tab << elapsedClock
               << tab << rate
               << tab << etaSeconds;
        file() << endl;

        // 2. Build progress bar
        string bar = "";
        const label filled = label((percent / 100.0) * barWidth_);
        for (label i = 0; i < barWidth_; ++i)
        {
            if (i < filled)
            {
                bar += "=";
            }
            else if (i == filled)
            {
                bar += ">";
            }
            else
            {
                bar += " ";
            }
        }

        // 3. Print formatted console output
        if (log_)
        {
            Info<< "Progress: [" << bar << "] "
                << std::fixed << std::setprecision(1) << percent << "% "
                << "(Time = " << curTime << " / " << endTime_ << "s) | "
                << "ETA: " << etaString << " | "
                << "Elapsed: " << formatTime(elapsedClock)
                << endl;
        }

        // 4. Optionally write machine-readable progress.json
        if (writeJson_)
        {
            fileName jsonPath = time_.path()/"postProcessing"/"progress.json";
            OFstream jsonFile(jsonPath);
            if (jsonFile.good())
            {
                jsonFile << "{\n"
                         << "  \"percent\": " << percent << ",\n"
                         << "  \"currentTime\": " << curTime << ",\n"
                         << "  \"startTime\": " << startTime_ << ",\n"
                         << "  \"endTime\": " << endTime_ << ",\n"
                         << "  \"elapsedSeconds\": " << elapsedClock << ",\n"
                         << "  \"etaSeconds\": " << etaSeconds << ",\n"
                         << "  \"rate\": " << rate << ",\n"
                         << "  \"completed\": " << (curTime >= endTime_ ? "true" : "false") << "\n"
                         << "}\n";
            }
        }
    }

    return true;
}


// ************************************************************************* //
