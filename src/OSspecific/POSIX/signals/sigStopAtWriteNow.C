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

\*---------------------------------------------------------------------------*/

#include "sigStopAtWriteNow.H"
#include "error.H"
#include "JobInfo.H"
#include "IOstreams.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
// Signal number to catch
int sigStopAtWriteNow::signal_
(
    debug::optimisationSwitch("stopAtWriteNowSignal", -1)
);
// Register re-reader
class addstopAtWriteNowSignalToOpt
:
    public ::Foam::simpleRegIOobject
{
public:
    addstopAtWriteNowSignalToOpt(const char* name)
    :
        ::Foam::simpleRegIOobject(Foam::debug::addOptimisationObject, name)
    {}
    virtual ~addstopAtWriteNowSignalToOpt()
    {}
    virtual void readData(Foam::Istream& is)
    {
        sigStopAtWriteNow::signal_ = readLabel(is);
        sigStopAtWriteNow::set(true);
    }
    virtual void writeData(Foam::Ostream& os) const
    {
        os << sigStopAtWriteNow::signal_;
    }
};
addstopAtWriteNowSignalToOpt addstopAtWriteNowSignalToOpt_
(
    "stopAtWriteNowSignal"
);
}


static Foam::Time const* runTimePtr_ = NULL;


struct sigaction Foam::sigStopAtWriteNow::oldAction_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sigStopAtWriteNow::sigHandler(int)
{
    // Reset old handling
    if (sigaction(signal_, &oldAction_, NULL) < 0)
    {
        FatalErrorIn
        (
            "Foam::sigStopAtWriteNow::sigHandler(int)"
        )   << "Cannot reset " << signal_ << " trapping"
            << abort(FatalError);
    }

    // Update jobInfo file
    jobInfo.signalEnd();

    Info<< "sigStopAtWriteNow :"
        << " setting up write and stop at end of the next iteration"
        << nl << endl;
    runTimePtr_->stopAt(Time::saWriteNow);

    //// Throw signal (to old handler)
    //raise(signal_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sigStopAtWriteNow::sigStopAtWriteNow(){}


Foam::sigStopAtWriteNow::sigStopAtWriteNow
(
    const bool verbose,
    const Time& runTime
)
{
    // Store runTime
    runTimePtr_ = &runTime;

    set(verbose);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sigStopAtWriteNow::~sigStopAtWriteNow()
{
    // Reset old handling
    if (signal_ > 0)
    {
        if (sigaction(signal_, &oldAction_, NULL) < 0)
        {
            FatalErrorIn
            (
                "Foam::sigStopAtWriteNow::~sigStopAtWriteNow()"
            )   << "Cannot reset " << signal_ << " trapping"
                << abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sigStopAtWriteNow::set(const bool verbose)
{
    if (signal_ > 0)
    {
        // Check that the signal is different from the writeNowSignal
        if (sigWriteNow::signal_ == signal_)
        {
            FatalErrorIn
            (
                "Foam::sigStopAtWriteNow::sigStopAtWriteNow"
                "(const bool, const Time&)"
            )   << "stopAtWriteNowSignal : " << signal_
                << " cannot be the same as the writeNowSignal."
                << " Please change this in the controlDict ("
                << findEtcFile("controlDict", false) << ")."
                << exit(FatalError);
        }


        struct sigaction newAction;
        newAction.sa_handler = sigHandler;
        newAction.sa_flags = SA_NODEFER;
        sigemptyset(&newAction.sa_mask);
        if (sigaction(signal_, &newAction, &oldAction_) < 0)
        {
            FatalErrorIn
            (
                "Foam::sigStopAtWriteNow::sigStopAtWriteNow"
                "(const bool, const Time&)"
            )   << "Cannot set " << signal_ << " trapping"
                << abort(FatalError);
        }

        if (verbose)
        {
            Info<< "sigStopAtWriteNow :"
                << " Enabling writing and stopping upon signal " << signal_
                << endl;
        }
    }
}


bool Foam::sigStopAtWriteNow::active() const
{
    return signal_ > 0;
}


// ************************************************************************* //
