/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "sigFpe.H"
#include "error.H"
#include "jobInfo.H"
#include "OSspecific.H"
#include "IOstreams.H"

#ifdef LINUX_GNUC
    #ifndef __USE_GNU
        #define __USE_GNU
    #endif
    #include <fenv.h>
    #include <malloc.h>
#endif

#include <limits>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

struct sigaction Foam::sigFpe::oldAction_;

void Foam::sigFpe::fillNan(UList<scalar>& lst)
{
    lst = std::numeric_limits<scalar>::signaling_NaN();
}

bool Foam::sigFpe::mallocNanActive_ = false;


#ifdef LINUX
extern "C"
{
    extern void* __libc_malloc(size_t size);

    // Override the GLIBC malloc to support mallocNan
    void* malloc(size_t size)
    {
        if (Foam::sigFpe::mallocNanActive_)
        {
            return Foam::sigFpe::mallocNan(size);
        }
        else
        {
            return __libc_malloc(size);
        }
    }
}

void* Foam::sigFpe::mallocNan(size_t size)
{
    // Call the low-level GLIBC malloc function
    void * result = __libc_malloc(size);

    // Initialize to signalling NaN
    UList<scalar> lst(reinterpret_cast<scalar*>(result), size/sizeof(scalar));
    sigFpe::fillNan(lst);

    return result;
}
#endif


#ifdef LINUX_GNUC
void Foam::sigFpe::sigHandler(int)
{
    // Reset old handling
    if (sigaction(SIGFPE, &oldAction_, nullptr) < 0)
    {
        FatalErrorInFunction
            << "Cannot reset SIGFPE trapping"
            << abort(FatalError);
    }

    // Update jobInfo file
    jobInfo_.signalEnd();

    error::printStack(Perr);

    // Throw signal (to old handler)
    raise(SIGFPE);
}
#endif


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sigFpe::sigFpe()
{
    oldAction_.sa_handler = nullptr;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sigFpe::~sigFpe()
{
    if (env("FOAM_SIGFPE"))
    {
        #ifdef LINUX_GNUC
        // Reset signal
        if
        (
            oldAction_.sa_handler
         && sigaction(SIGFPE, &oldAction_, nullptr) < 0
        )
        {
            FatalErrorInFunction
                << "Cannot reset SIGFPE trapping"
                << abort(FatalError);
        }
        #endif
    }

    if (env("FOAM_SETNAN"))
    {
        #ifdef LINUX
        // Disable initialization to NaN
        mallocNanActive_ = false;
        #endif
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sigFpe::set(const bool verbose)
{
    if (oldAction_.sa_handler)
    {
        FatalErrorInFunction
            << "Cannot call sigFpe::set() more than once"
            << abort(FatalError);
    }

    if (env("FOAM_SIGFPE"))
    {
        bool supported = false;

        #ifdef LINUX_GNUC
        supported = true;

        feenableexcept
        (
            FE_DIVBYZERO
          | FE_INVALID
          | FE_OVERFLOW
        );

        struct sigaction newAction;
        newAction.sa_handler = sigHandler;
        newAction.sa_flags = SA_NODEFER;
        sigemptyset(&newAction.sa_mask);
        if (sigaction(SIGFPE, &newAction, &oldAction_) < 0)
        {
            FatalErrorInFunction
                << "Cannot set SIGFPE trapping"
                << abort(FatalError);
        }
        #endif

        if (verbose)
        {
            if (supported)
            {
                Info<< "sigFpe : Enabling floating point exception trapping"
                    << " (FOAM_SIGFPE)." << endl;
            }
            else
            {
                Info<< "sigFpe : Floating point exception trapping"
                    << " - not supported on this platform" << endl;
            }
        }
    }


    if (env("FOAM_SETNAN"))
    {
        #ifdef LINUX
        mallocNanActive_ = true;
        #endif

        if (verbose)
        {
            if (mallocNanActive_)
            {
                Info<< "SetNaN : Initialising allocated memory to NaN"
                    << " (FOAM_SETNAN)." << endl;
            }
            else
            {
                Info<< "SetNaN : Initialise allocated memory to NaN"
                    << " - not supported on this platform" << endl;
            }
        }
    }
}


// ************************************************************************* //
