/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "fileMonitor.H"
#include "IOstreams.H"
#include "Pstream.H"
#include "PackedList.H"
#include "PstreamReduceOps.H"
#include "OSspecific.H"
#include "regIOobject.H"     // for fileModificationSkew symbol

#ifdef FOAM_USE_INOTIFY
#   include <unistd.h>
#   include <sys/inotify.h>
#   include <sys/ioctl.h>
#   include <errno.h>
#   define EVENT_SIZE  ( sizeof (struct inotify_event) )
#   define EVENT_LEN   (EVENT_SIZE + 16)
#   define EVENT_BUF_LEN     ( 1024 * EVENT_LEN )
#else
#   include "OSspecific.H"
#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::NamedEnum<Foam::fileMonitor::fileState, 3>
    Foam::fileMonitor::fileStateNames_;

namespace Foam
{
    defineTypeNameAndDebug(fileMonitor, 0);

    template<>
    const char* Foam::NamedEnum
    <
        Foam::fileMonitor::fileState,
        3
    >::names[] =
    {
        "unmodified",
        "modified",
        "deleted"
    };

    //- Reduction operator for PackedList of fileState
    class reduceFileStates
    {
        public:
        unsigned int operator()(const unsigned int x, const unsigned int y)
        const
        {
            // x,y are sets of 2bits representing fileState

            unsigned int mask = 3u;
            unsigned int shift = 0;
            unsigned int result = 0;

            while (mask)
            {
                // Combine state
                unsigned int xState = (x & mask) >> shift;
                unsigned int yState = (y & mask) >> shift;

                // Combine and add to result. Combine is such that UNMODIFIED
                // wins.
                unsigned int state = min(xState, yState);
                result |= (state << shift);

                shift += 2;
                mask <<= 2;
            }
            return result;
        }
    };

    //- Combine operator for PackedList of fileState
    class combineReduceFileStates
    {
        public:
        void operator()(unsigned int& x, const unsigned int y) const
        {
            x = reduceFileStates()(x, y);
        }
    };



    //-  Internal tracking via stat(3p) or inotify(7)
    class fileMonitorWatcher
    {
    public:

        const bool useInotify_;

        // For inotify

            //- File descriptor for the inotify instance
            int inotifyFd_;

            //- Current watchIDs and corresponding directory id
            DynamicList<label> dirWatches_;
            DynamicList<fileName> dirFiles_;

        // For stat

            //- From watch descriptor to modified time
            DynamicList<time_t> lastMod_;



        //- initialise inotify
        inline fileMonitorWatcher(const bool useInotify, const label sz = 20)
        :
            useInotify_(useInotify),
            inotifyFd_(-1)
        {
            if (useInotify_)
            {
                #ifdef FOAM_USE_INOTIFY
                inotifyFd_ = inotify_init();
                dirWatches_.setCapacity(sz);
                dirFiles_.setCapacity(sz);

                if (inotifyFd_ < 0)
                {
                    static bool hasWarned = false;
                    if (!hasWarned)
                    {
                        hasWarned = true;
                        WarningIn("fileMonitorWatcher(const bool, const label)")
                            << "Failed allocating an inotify descriptor : "
                            << string(strerror(errno)) << endl
                            << "    Please increase the number of allowable "
                            << "inotify instances" << endl
                            << "    (/proc/sys/fs/inotify/max_user_instances"
                            << " on Linux)" << endl
                            << "    , switch off runTimeModifiable." << endl
                            << "    or compile this file without "
                            << "FOAM_USE_INOTIFY"
                            << " to use time stamps instead of inotify." << endl
                            << "    Continuing without additional file"
                            << " monitoring."
                            << endl;
                    }
                }
                #else
                    FatalErrorIn("fileMonitorWatcher(const bool, const label)")
                        << "You selected inotify but this file was compiled"
                        << " without FOAM_USE_INOTIFY"
                        << " Please select another fileModification test method"
                        << exit(FatalError);
                #endif
            }
            else
            {
                lastMod_.setCapacity(sz);
            }
        }

        //- remove all watches
        inline ~fileMonitorWatcher()
        {
            #ifdef FOAM_USE_INOTIFY
            if (useInotify_ && inotifyFd_ >= 0)
            {
                forAll(dirWatches_, i)
                {
                    if (dirWatches_[i] >= 0)
                    {
                        if (inotify_rm_watch(inotifyFd_, int(dirWatches_[i])))
                        {
                            WarningIn("fileMonitor::~fileMonitor()")
                                << "Failed deleting directory watch "
                                << dirWatches_[i] << endl;
                        }
                    }
                }
            }
            #endif
        }

        inline bool addWatch(const label watchFd, const fileName& fName)
        {
            if (useInotify_)
            {
                if (inotifyFd_ < 0)
                {
                    return false;
                }

                #ifdef FOAM_USE_INOTIFY
                // Add/retrieve watch on directory containing file.
                // Note that fName might be non-existing in special situations
                // (master-only reading for IODictionaries)

                const fileName dir = fName.path();

                label dirWatchID = -1;
                if (isDir(dir))
                {
                    dirWatchID = inotify_add_watch
                    (
                        inotifyFd_,
                        dir.c_str(),
                        IN_CLOSE_WRITE
                    );

                    if (dirWatchID < 0)
                    {
                        FatalErrorIn("addWatch(const label, const fileName&)")
                            << "Failed adding watch " << watchFd
                            << " to directory " << fName << " due to "
                            << string(strerror(errno))
                            << exit(FatalError);
                    }
                }

                if (watchFd < dirWatches_.size() && dirWatches_[watchFd] != -1)
                {
                    // Reuse of watchFd : should have dir watchID set to -1.
                    FatalErrorIn("addWatch(const label, const fileName&)")
                        << "Problem adding watch " << watchFd
                        << " to file " << fName
                        << abort(FatalError);
                }

                dirWatches_(watchFd) = dirWatchID;
                dirFiles_(watchFd) = fName.name();
                #endif
            }
            else
            {
                if (watchFd < lastMod_.size() && lastMod_[watchFd] != 0)
                {
                    // Reuse of watchFd : should have lastMod set to 0.
                    FatalErrorIn("addWatch(const label, const fileName&)")
                        << "Problem adding watch " << watchFd
                        << " to file " << fName
                        << abort(FatalError);
                }

                lastMod_(watchFd) = lastModified(fName);
            }

            return true;
        }

        inline bool removeWatch(const label watchFd)
        {
            if (useInotify_)
            {
                if (inotifyFd_ < 0)
                {
                    return false;
                }

                dirWatches_[watchFd] = -1;
            }
            else
            {
                lastMod_[watchFd] = 0;
            }
            return true;
        }

    };
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fileMonitor::checkFiles() const
{
    if (useInotify_)
    {
        #ifdef FOAM_USE_INOTIFY
        // Large buffer for lots of events
        char buffer[EVENT_BUF_LEN];

        while (true)
        {
            struct timeval zeroTimeout = {0, 0};

            //- Pre-allocated structure containing file descriptors
            fd_set fdSet;
            // Add notify descriptor to select fd_set
            FD_ZERO(&fdSet);
            FD_SET(watcher_->inotifyFd_, &fdSet);

            int ready = select
            (
                watcher_->inotifyFd_+1,     // num filedescriptors in fdSet
                &fdSet,             // fdSet with only inotifyFd
                NULL,               // No writefds
                NULL,               // No errorfds
                &zeroTimeout        // eNo timeout
            );

            if (ready < 0)
            {
                FatalErrorIn("fileMonitor::checkFiles()")
                    << "Problem in issuing select."
                    << abort(FatalError);
            }
            else if (FD_ISSET(watcher_->inotifyFd_, &fdSet))
            {
                // Read events
                ssize_t nBytes = read
                (
                    watcher_->inotifyFd_,
                    buffer,
                    EVENT_BUF_LEN
                );

                if (nBytes < 0)
                {
                    FatalErrorIn("fileMonitor::checkFiles()")
                        << "read of " << watcher_->inotifyFd_
                        << " failed with " << label(nBytes)
                        << abort(FatalError);
                }

                // Go through buffer, consuming events
                int i = 0;
                while (i < nBytes)
                {
                    const struct inotify_event* inotifyEvent =
                        reinterpret_cast<const struct inotify_event*>
                        (
                            &buffer[i]
                        );

                    //Pout<< "watchFd:" << inotifyEvent->wd << nl
                    //    << "mask:" << inotifyEvent->mask << nl
                    //  << endl;
                    //Pout<< "file:" << fileName(inotifyEvent->name) << endl;
                    //Pout<< "len:" << inotifyEvent->len << endl;

                    if
                    (
                        (inotifyEvent->mask & IN_CLOSE_WRITE)
                     && inotifyEvent->len
                    )
                    {
                        // Search for file
                        forAll(watcher_->dirWatches_, i)
                        {
                            label id = watcher_->dirWatches_[i];
                            if
                            (
                                id == inotifyEvent->wd
                             && inotifyEvent->name == watcher_->dirFiles_[i]
                            )
                            {
                                // Correct directory and name
                                localState_[i] = MODIFIED;
                            }
                        }
                    }

                    i += EVENT_SIZE + inotifyEvent->len;
                }
            }
            else
            {
                // No data
                return;
            }
        }
        #endif
    }
    else
    {
        forAll(watcher_->lastMod_, watchFd)
        {
            time_t oldTime = watcher_->lastMod_[watchFd];

            if (oldTime != 0)
            {
                const fileName& fName = watchFile_[watchFd];
                time_t newTime = lastModified(fName);

                if (newTime == 0)
                {
                    localState_[watchFd] = DELETED;
                }
                else
                {
                    if (newTime > (oldTime + regIOobject::fileModificationSkew))
                    {
                        localState_[watchFd] = MODIFIED;
                    }
                    else
                    {
                        localState_[watchFd] = UNMODIFIED;
                    }
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::fileMonitor::fileMonitor(const bool useInotify)
:
    useInotify_(useInotify),
    localState_(20),
    state_(20),
    watchFile_(20),
    freeWatchFds_(2),
    watcher_(new fileMonitorWatcher(useInotify_, 20))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fileMonitor::~fileMonitor()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Note: fName might not exist (on slaves if in master-only mode for
// regIOobject)
Foam::label Foam::fileMonitor::addWatch(const fileName& fName)
{
    label watchFd;

    label sz = freeWatchFds_.size();

    if (sz)
    {
        watchFd = freeWatchFds_[sz-1];
        freeWatchFds_.setSize(sz-1);
    }
    else
    {
        watchFd = state_.size();
    }

    watcher_->addWatch(watchFd, fName);

    if (debug)
    {
        Pout<< "fileMonitor : added watch " << watchFd << " on file "
            << fName << endl;
    }

    if (watchFd < 0)
    {
        WarningIn("fileMonitor::addWatch(const fileName&)")
            << "could not add watch for file " << fName << endl;
    }
    else
    {
        localState_(watchFd) = UNMODIFIED;
        state_(watchFd) = UNMODIFIED;
        watchFile_(watchFd) = fName;
    }
    return watchFd;
}


bool Foam::fileMonitor::removeWatch(const label watchFd)
{
    if (debug)
    {
        Pout<< "fileMonitor : removing watch " << watchFd << " on file "
            << watchFile_[watchFd] << endl;
    }

    freeWatchFds_.append(watchFd);
    return watcher_->removeWatch(watchFd);
}


const Foam::fileName& Foam::fileMonitor::getFile(const label watchFd) const
{
    return watchFile_[watchFd];
}


Foam::fileMonitor::fileState Foam::fileMonitor::getState(const label watchFd)
const
{
    return state_[watchFd];
}


void Foam::fileMonitor::updateStates
(
    const bool masterOnly,
    const bool syncPar
) const
{
    if (Pstream::master() || !masterOnly)
    {
        // Update the localState_
        checkFiles();
    }

    if (syncPar)
    {
        // Pack local state (might be on master only)
        PackedList<2> stats(state_.size(), MODIFIED);
        if (Pstream::master() || !masterOnly)
        {
            forAll(state_, watchFd)
            {
                stats[watchFd] = static_cast<unsigned int>
                (
                    localState_[watchFd]
                );
            }
        }


        // Scatter or reduce to synchronise state
        if (masterOnly)
        {
            // Scatter
            if (stats.storage().size() == 1)
            {
                Pstream::scatter(stats.storage()[0]);
            }
            else
            {
                Pstream::listCombineScatter(stats.storage());
            }
        }
        else
        {
            // Reduce
            if (stats.storage().size() == 1)
            {
                // Optimisation valid for most cases.
                reduce(stats.storage()[0], reduceFileStates());
            }
            else
            {
                Pstream::listCombineGather
                (
                    stats.storage(),
                    combineReduceFileStates()
                );
            }
        }


        // Update synchronised state
        forAll(state_, watchFd)
        {
            // Assign synchronised state
            unsigned int stat = stats[watchFd];
            state_[watchFd] = fileState(stat);

            if (!masterOnly)
            {
                // Give warning for inconsistent state
                if (state_[watchFd] != localState_[watchFd])
                {
                    if (debug)
                    {
                        Pout<< "fileMonitor : Delaying reading "
                            << watchFile_[watchFd]
                            << " due to inconsistent "
                               "file time-stamps between processors"
                            << endl;
                    }

                    WarningIn
                    (
                        "fileMonitor::updateStates"
                        "(const bool, const bool) const"
                    )   << "Delaying reading " << watchFile_[watchFd]
                        << " due to inconsistent "
                           "file time-stamps between processors" << endl;
                }
            }
        }
    }
    else
    {
        state_ = localState_;
    }
}


void Foam::fileMonitor::setUnmodified(const label watchFd)
{
    state_[watchFd] = UNMODIFIED;
    localState_[watchFd] = UNMODIFIED;

    if (!useInotify_)
    {
        watcher_->lastMod_[watchFd] = lastModified(watchFile_[watchFd]);
    }
}


// ************************************************************************* //
