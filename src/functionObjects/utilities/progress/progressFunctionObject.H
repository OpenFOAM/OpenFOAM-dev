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

Class
    Foam::functionObjects::progress

Description
    Calculates and reports real-time simulation progress, percentage completion,
    solving speed, and estimated time of arrival (ETA).

    Example of function object specification:
    \verbatim
    progress
    {
        type            progress;
        libs            ("libutilityFunctionObjects.so");
        writeControl    timeStep;
        writeInterval   1;
        log             true;
        writeJson       true;
    }
    \endverbatim

    Output data is written to the dir postProcessing/progress/\<timeDir\>/
    and optionally to standard log output and a machine-readable JSON file.

See also
    Foam::functionObject
    Foam::functionObjects::logFiles

SourceFiles
    progressFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef progressFunctionObject_H
#define progressFunctionObject_H

#include "functionObject.H"
#include "logFiles.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class progress Declaration
\*---------------------------------------------------------------------------*/

class progress
:
    public functionObject,
    public logFiles
{
    // Private Data

        //- Simulation start time
        scalar startTime_;

        //- Simulation end time
        scalar endTime_;

        //- Whether to print progress to screen/log
        bool log_;

        //- Whether to write progress.json
        bool writeJson_;

        //- Width of progress bar in characters
        label barWidth_;


    // Private Member Functions

        //- Format seconds into a human-readable string (hh:mm:ss or mm:ss)
        string formatTime(const scalar seconds) const;

        //- Output file header information
        virtual void writeFileHeader(const label i);


public:

    //- Runtime type information
    TypeName("progress");


    // Constructors

        //- Construct from Time and dictionary
        progress
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        progress(const progress&) = delete;


    //- Destructor
    virtual ~progress();


    // Member Functions

        //- Read the controls
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const
        {
            return wordList::null();
        }

        //- Execute function object
        virtual bool execute();

        //- Write progress information
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const progress&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
