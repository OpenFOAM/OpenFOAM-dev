/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

Class
    Foam::sigFpe

Description
    Set up trapping for floating point exceptions (signal FPE).

    Controlled by two env vars:
      - \par FOAM_SIGFPE
        Exception trapping
      - \par FOAM_SETNAN
        Initialisation of all malloced memory to NaN. If FOAM_SIGFPE
        also set, this will cause usage of uninitialised scalars to trigger
        an abort.

SourceFiles
    sigFpe.C

\*---------------------------------------------------------------------------*/

#ifndef sigFpe_H
#define sigFpe_H

#include <signal.h>

#if defined(linux) || defined(linux64) || \
    defined(linuxARM7) || defined(linuxPPC64) || defined(linuxPPC64le)
    #define LINUX
#endif

#if defined(LINUX) && defined(__GNUC__)
    #define LINUX_GNUC
#endif

#include "UList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class sigFpe Declaration
\*---------------------------------------------------------------------------*/

class sigFpe
{
    // Private Data

        //- Saved old signal trapping setting
        static struct sigaction oldAction_;


    // Static Data Members

        #ifdef LINUX_GNUC
        //- Handler for caught signals
        static void sigHandler(int);
        #endif


public:


    // Constructors

        //- Construct null
        sigFpe();


    //- Destructor
    ~sigFpe();


    // Member Functions

        //- Activate SIGFPE signal handler when FOAM_SIGFPE is %set
        //  Fill memory with NaN when FOAM_SETNAN is %set
        void set(const bool verbose);

        //- Flag to indicate mallocNan is enabled
        static bool mallocNanActive_;

        #ifdef LINUX
        //- Malloc function which initialises to NaN
        static void* mallocNan(size_t size);
        #endif

        //- Fill block of data with NaN
        static void fillNan(UList<scalar>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
