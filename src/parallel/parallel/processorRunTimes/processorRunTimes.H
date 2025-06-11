/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2025 OpenFOAM Foundation
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
    Foam::processorRunTimes

Description

SourceFiles
    processorRunTimes.C

\*---------------------------------------------------------------------------*/

#ifndef processorRunTimes_H
#define processorRunTimes_H

#include "Time.H"
#include "PtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class processorRunTimes Declaration
\*---------------------------------------------------------------------------*/

class processorRunTimes
{
public:

    // Public Enumerations

        //- Enumeration for how the number of processors is determined
        enum class nProcsFrom
        {
            decomposeParDict,
            fileHandler
        };


private:

    // Private Data

        //- The complete run time
        Time completeRunTime_;

        //- Processor run times
        PtrList<Time> procRunTimes_;

        //- Processor zero time that can be used/queried when there is no
        //  decomposition. Only allocated when procRunTimes_ is empty.
        autoPtr<Time> proc0RunTime_;


    // Private Member Functions

        //- Construct the processor run times
        static PtrList<Time> procRunTimes
        (
            const Time& completeRunTime,
            const nProcsFrom
        );

        //- Construct the processor0 run time (if necessary)
        static autoPtr<Time> proc0RunTime
        (
            const Time& completeRunTime,
            const PtrList<Time>& procRunTimes
        );

        //- Access the processor0 run time
        inline Time& proc0TimeRef()
        {
            return proc0RunTime_.valid() ? proc0RunTime_() : procRunTimes_[0];
        }


public:

    // Constructors

        //- Construct from name and arguments
        processorRunTimes
        (
            const word& name,
            const argList& args,
            const bool enableFunctionObjects = true,
            const nProcsFrom = nProcsFrom::decomposeParDict
        );

        //- Construct from name and directories
        processorRunTimes
        (
            const word& name,
            const fileName& rootPath,
            const fileName& caseName,
            const bool enableFunctionObjects = true,
            const nProcsFrom = nProcsFrom::decomposeParDict
        );


    //- Destructor
    ~processorRunTimes();


    // Member Functions

        //- Access the complete run time
        inline const Time& completeTime() const
        {
            return completeRunTime_;
        }

        //- Access the processor run times
        inline const PtrList<Time>& procTimes() const
        {
            return procRunTimes_;
        }

        //- Return the number of processors
        inline label nProcs() const
        {
            return procRunTimes_.size();
        }

        //- Access the processor0 run time
        inline const Time& proc0Time() const
        {
            return proc0RunTime_.valid() ? proc0RunTime_() : procRunTimes_[0];
        }

        //- Set the time
        void setTime(const instant& inst, const label newIndex);

        //- Select the time
        instantList selectComplete(const argList& args);

        //- Select the time
        instantList selectProc(const argList& args);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
