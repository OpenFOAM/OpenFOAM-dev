/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "outputFilterControl.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* NamedEnum<outputFilterControl::writeControls, 8>::
    names[] =
    {
        "timeStep",
        "writeTime",
        "outputTime",
        "adjustableRunTime",
        "runTime",
        "clockTime",
        "cpuTime",
        "none"
    };
}

const Foam::NamedEnum<Foam::outputFilterControl::writeControls, 8>
    Foam::outputFilterControl::writeControlNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::outputFilterControl::outputFilterControl
(
    const Time& t,
    const dictionary& dict,
    const word& prefix
)
:
    time_(t),
    prefix_(prefix),
    writeControl_(ocTimeStep),
    outputInterval_(0),
    writeInterval_(-1),
    outputTimeLastDump_(0)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::outputFilterControl::~outputFilterControl()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::outputFilterControl::read(const dictionary& dict)
{
    word controlName(prefix_ + "Control");
    word intervalName(prefix_ + "Interval");

    // For backward compatibility support the deprecated 'outputControl' option
    // now superseded by 'writeControl' for compatibility with Time
    if (prefix_ == "write" && dict.found("outputControl"))
    {
        IOWarningInFunction(dict)
            << "Using deprecated 'outputControl'" << nl
            << "    Please use 'writeControl' with 'writeInterval'"
            << endl;

        // Change to the old names for this option
        controlName = "outputControl";
        intervalName = "outputInterval";
    }

    if (dict.found(controlName))
    {
        writeControl_ = writeControlNames_.read(dict.lookup(controlName));
    }
    else
    {
        writeControl_ = ocTimeStep;
    }

    switch (writeControl_)
    {
        case ocTimeStep:
        {
            outputInterval_ = dict.lookupOrDefault<label>(intervalName, 0);
            break;
        }

        case ocWriteTime:
        case ocOutputTime:
        {
            outputInterval_ = dict.lookupOrDefault<label>(intervalName, 1);
            break;
        }

        case ocClockTime:
        case ocRunTime:
        case ocCpuTime:
        case ocAdjustableRunTime:
        {
            writeInterval_ = readScalar(dict.lookup(intervalName));
            break;
        }

        default:
        {
            break;
        }
    }
}


bool Foam::outputFilterControl::output()
{
    switch (writeControl_)
    {
        case ocTimeStep:
        {
            return
            (
                (outputInterval_ <= 1)
             || !(time_.timeIndex() % outputInterval_)
            );
            break;
        }

        case ocWriteTime:
        case ocOutputTime:
        {
            if (time_.outputTime())
            {
                outputTimeLastDump_++;
                return !(outputTimeLastDump_ % outputInterval_);
            }
            break;
        }

        case ocRunTime:
        case ocAdjustableRunTime:
        {
            label outputIndex = label
            (
                (
                    (time_.value() - time_.startTime().value())
                  + 0.5*time_.deltaTValue()
                )
               /writeInterval_
            );

            if (outputIndex > outputTimeLastDump_)
            {
                outputTimeLastDump_ = outputIndex;
                return true;
            }
            break;
        }

        case ocCpuTime:
        {
            label outputIndex = label
            (
                returnReduce(time_.elapsedCpuTime(), maxOp<double>())
               /writeInterval_
            );
            if (outputIndex > outputTimeLastDump_)
            {
                outputTimeLastDump_ = outputIndex;
                return true;
            }
            break;
        }

        case ocClockTime:
        {
            label outputIndex = label
            (
                returnReduce(label(time_.elapsedClockTime()), maxOp<label>())
               /writeInterval_
            );
            if (outputIndex > outputTimeLastDump_)
            {
                outputTimeLastDump_ = outputIndex;
                return true;
            }
            break;
        }

        case ocNone:
        {
            return false;
        }

        default:
        {
            FatalErrorInFunction
                << "Undefined output control: "
                << writeControlNames_[writeControl_] << nl
                << abort(FatalError);
            break;
        }
    }

    return false;
}


// ************************************************************************* //
