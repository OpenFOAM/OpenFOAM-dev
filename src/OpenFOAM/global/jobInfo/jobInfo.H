/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::jobInfo

Description
    Helper class for recording information about run/finished jobs.

    Writes the following files:
      - $FOAM_JOB_DIR/runningJobs
      - $FOAM_JOB_DIR/finishedJobs

SourceFiles
    jobInfo.C

\*---------------------------------------------------------------------------*/

#ifndef jobInfo_H
#define jobInfo_H

#include "dictionary.H"
#include "fileName.H"
#include "cpuTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class jobInfo Declaration
\*---------------------------------------------------------------------------*/

class jobInfo
:
    public dictionary
{
    // Private Data

        mutable fileName runningJobPath_;
        mutable fileName finishedJobPath_;
        cpuTime cpuTime_;

    // Private Member Functions

        bool write(Ostream&) const;
        void end(const word& terminationType);


public:

    static bool constructed;
    static bool writeJobControl;
    static bool writeJobInfo;

    // Constructors

        //- Construct null
        jobInfo();


    //- Destructor
    ~jobInfo();


    // Member Functions

        void write(const word& executable, const fileName& casePath) const;

        void end();

        void exit();

        void abort();

        void signalEnd() const;
};


extern jobInfo jobInfo_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
