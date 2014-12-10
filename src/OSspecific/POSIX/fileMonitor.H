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

Class
    fileMonitor

Description
    Checking for changes to files.

Note
    The default is to use stat to get the timestamp.

    Compile with FOAM_USE_INOTIFY to use the inotify
    (Linux specific, since 2.6.13) framework. The problem is that inotify does
    not work on nfs3 mounted directories!!

SourceFiles
    fileMonitor.C

\*---------------------------------------------------------------------------*/

#ifndef fileMonitor_H
#define fileMonitor_H

#include <sys/types.h>
#include "NamedEnum.H"
#include "className.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fileMonitor;
class fileMonitorWatcher;

/*---------------------------------------------------------------------------*\
                         Class fileMonitor Declaration
\*---------------------------------------------------------------------------*/

class fileMonitor
{

public:

    // Public data types

        //- Enumeration defining the file state.
        enum fileState
        {
            UNMODIFIED = 0,
            MODIFIED = 1,
            DELETED = 2
        };

        static const NamedEnum<fileState, 3> fileStateNames_;

private:
    // Private data

        //- Whether to use inotify (requires -DFOAM_USE_INOTIFY, see above)
        const bool useInotify_;

        //- State for all watchFds based on local files
        mutable DynamicList<fileState> localState_;

        //- State for all watchFds - synchronised
        mutable DynamicList<fileState> state_;

        //- Filename for all watchFds
        DynamicList<fileName> watchFile_;

        //- Free watchFds
        DynamicList<label> freeWatchFds_;

        //- Watch mechanism (stat or inotify)
        mutable autoPtr<fileMonitorWatcher> watcher_;


    // Private Member Functions

        //- Update localState_ from any events.
        void checkFiles() const;

        //- Disallow default bitwise copy construct
        fileMonitor(const fileMonitor&);

        //- Disallow default bitwise assignment
        void operator=(const fileMonitor&);


public:

    // Declare name of the class and its debug switch
    ClassName("fileMonitor");

    // Constructors

        //- Construct null
        fileMonitor(const bool useInotify);


    //- Destructor
    ~fileMonitor();


    // Member Functions

        //- Add file to watch. Return watch descriptor
        label addWatch(const fileName&);

        //- Remove file to watch. Return true if successful
        bool removeWatch(const label watchFd);

        //- Get name of file being watched
        const fileName& getFile(const label watchFd) const;

        //- Check state using handle
        fileState getState(const label watchFd) const;

        //- Check state of all files. Updates state_.
        void updateStates
        (
            const bool masterOnly,
            const bool syncPar
        ) const;

        //- Reset state (e.g. after having read it) using handle
        void setUnmodified(const label watchFd);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
