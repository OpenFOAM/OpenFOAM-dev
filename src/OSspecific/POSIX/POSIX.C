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

Description
    POSIX versions of the functions declared in OSspecific.H

\*---------------------------------------------------------------------------*/

#ifdef solarisGcc
    #define _SYS_VNODE_H
#endif

#include "OSspecific.H"
#include "POSIX.H"
#include "foamVersion.H"
#include "fileName.H"
#include "fileStat.H"
#include "timer.H"
#include "IFstream.H"
#include "DynamicList.H"
#include "HashSet.H"
#include "IOstreams.H"
#include "Pstream.H"

#include <fstream>
#include <cstdlib>
#include <cctype>

#include <stdio.h>
#include <unistd.h>
#include <dirent.h>
#include <pwd.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/socket.h>
#include <netdb.h>
#include <dlfcn.h>
#include <link.h>

#include <netinet/in.h>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(POSIX, 0);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

pid_t Foam::pid()
{
    return ::getpid();
}


pid_t Foam::ppid()
{
    return ::getppid();
}


pid_t Foam::pgid()
{
    return ::getpgrp();
}


bool Foam::env(const word& envName)
{
    return ::getenv(envName.c_str()) != nullptr;
}


Foam::string Foam::getEnv(const word& envName)
{
    char* env = ::getenv(envName.c_str());

    if (env)
    {
        return string(env);
    }
    else
    {
        // Return null-constructed string rather than string::null
        // to avoid cyclic dependencies in the construction of globals
        return string();
    }
}


bool Foam::setEnv
(
    const word& envName,
    const std::string& value,
    const bool overwrite
)
{
    return setenv(envName.c_str(), value.c_str(), overwrite) == 0;
}


Foam::string Foam::hostName(bool full)
{
    char buf[128];
    ::gethostname(buf, sizeof(buf));

    // Implementation as per hostname from net-tools
    if (full)
    {
        struct hostent *hp = ::gethostbyname(buf);
        if (hp)
        {
            return hp->h_name;
        }
    }

    return buf;
}


Foam::string Foam::domainName()
{
    char buf[128];
    ::gethostname(buf, sizeof(buf));

    // Implementation as per hostname from net-tools
    struct hostent *hp = ::gethostbyname(buf);
    if (hp)
    {
        char *p = ::strchr(hp->h_name, '.');
        if (p)
        {
            ++p;
            return p;
        }
    }

    return string::null;
}


Foam::string Foam::userName()
{
    struct passwd* pw = ::getpwuid(::getuid());

    if (pw != nullptr)
    {
        return pw->pw_name;
    }
    else
    {
        return string::null;
    }
}


bool Foam::isAdministrator()
{
    return (::geteuid() == 0);
}


Foam::fileName Foam::home()
{
    char* env = ::getenv("HOME");

    if (env != nullptr)
    {
        return fileName(env);
    }
    else
    {
        struct passwd* pw = ::getpwuid(getuid());

        if (pw != nullptr)
        {
            return pw->pw_dir;
        }
        else
        {
            return fileName::null;
        }
    }
}


Foam::fileName Foam::home(const string& userName)
{
    struct passwd* pw;

    if (userName.size())
    {
        pw = ::getpwnam(userName.c_str());
    }
    else
    {
        char* env = ::getenv("HOME");

        if (env != nullptr)
        {
            return fileName(env);
        }

        pw = ::getpwuid(::getuid());
    }

    if (pw != nullptr)
    {
        return pw->pw_dir;
    }
    else
    {
        return fileName::null;
    }
}


Foam::fileName Foam::cwd()
{
    label pathLengthLimit = POSIX::pathLengthChunk;
    List<char> path(pathLengthLimit);

    // Resize path if getcwd fails with an ERANGE error
    while(pathLengthLimit == path.size())
    {
        if (::getcwd(path.data(), path.size()))
        {
            return path.data();
        }
        else if(errno == ERANGE)
        {
            // Increment path length up to the pathLengthMax limit
            if
            (
                (pathLengthLimit += POSIX::pathLengthChunk)
             >= POSIX::pathLengthMax
            )
            {
                FatalErrorInFunction
                    << "Attempt to increase path length beyond limit of "
                    << POSIX::pathLengthMax
                    << exit(FatalError);
            }

            path.setSize(pathLengthLimit);
        }
        else
        {
            break;
        }
    }

    FatalErrorInFunction
        << "Couldn't get the current working directory"
        << exit(FatalError);

    return fileName::null;
}


bool Foam::chDir(const fileName& dir)
{
    return ::chdir(dir.c_str()) == 0;
}


bool Foam::mkDir(const fileName& filePath, mode_t mode)
{
    if (POSIX::debug)
    {
        Pout<< FUNCTION_NAME << " : filePath:" << filePath << " mode:" << mode
            << endl;
        if ((POSIX::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    // Empty names are meaningless
    if (filePath.empty())
    {
        return false;
    }

    // Construct instance path directory if does not exist
    if (::mkdir(filePath.c_str(), mode) == 0)
    {
        // Directory made OK so return true
        return true;
    }
    else
    {
        switch (errno)
        {
            case EPERM:
            {
                FatalErrorInFunction
                    << "The filesystem containing " << filePath
                    << " does not support the creation of directories."
                    << exit(FatalError);

                return false;
            }

            case EEXIST:
            {
                // Directory already exists so simply return true
                return true;
            }

            case EFAULT:
            {
                FatalErrorInFunction
                    << "" << filePath
                    << " points outside your accessible address space."
                    << exit(FatalError);

                return false;
            }

            case EACCES:
            {
                FatalErrorInFunction
                    << "The parent directory does not allow write "
                       "permission to the process,"<< nl
                    << "or one of the directories in " << filePath
                    << " did not allow search (execute) permission."
                    << exit(FatalError);

                return false;
            }

            case ENAMETOOLONG:
            {
                FatalErrorInFunction
                    << "" << filePath << " is too long."
                    << exit(FatalError);

                return false;
            }

            case ENOENT:
            {
                // Part of the path does not exist so try to create it
                if (filePath.path().size() && mkDir(filePath.path(), mode))
                {
                    return mkDir(filePath, mode);
                }
                else
                {
                    FatalErrorInFunction
                        << "Couldn't create directory " << filePath
                        << exit(FatalError);

                    return false;
                }
            }

            case ENOTDIR:
            {
                FatalErrorInFunction
                    << "A component used as a directory in " << filePath
                    << " is not, in fact, a directory."
                    << exit(FatalError);

                return false;
            }

            case ENOMEM:
            {
                FatalErrorInFunction
                    << "Insufficient kernel memory was available to make "
                       "directory " << filePath << '.'
                    << exit(FatalError);

                return false;
            }

            case EROFS:
            {
                FatalErrorInFunction
                    << "" << filePath
                    << " refers to a file on a read-only filesystem."
                    << exit(FatalError);

                return false;
            }

            case ELOOP:
            {
                FatalErrorInFunction
                    << "Too many symbolic links were encountered in resolving "
                    << filePath << '.'
                    << exit(FatalError);

                return false;
            }

            case ENOSPC:
            {
                FatalErrorInFunction
                    << "The device containing " << filePath
                    << " has no room for the new directory or "
                    << "the user's disk quota is exhausted."
                    << exit(FatalError);

                return false;
            }

            default:
            {
                FatalErrorInFunction
                    << "Couldn't create directory " << filePath
                    << exit(FatalError);

                return false;
            }
        }
    }
}


bool Foam::chMod(const fileName& name, const mode_t m)
{
    if (POSIX::debug)
    {
        Pout<< FUNCTION_NAME << " : name:" << name << endl;
        if ((POSIX::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }
    return ::chmod(name.c_str(), m) == 0;
}


mode_t Foam::mode
(
    const fileName& name,
    const bool checkVariants,
    const bool followLink
)
{
    if (POSIX::debug)
    {
        Pout<< FUNCTION_NAME << " : name:" << name << endl;
        if ((POSIX::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }
    fileStat fileStatus(name, checkVariants, followLink);
    if (fileStatus.isValid())
    {
        return fileStatus.status().st_mode;
    }
    else
    {
        return 0;
    }
}


Foam::fileType Foam::type
(
    const fileName& name,
    const bool checkVariants,
    const bool followLink
)
{
    if (POSIX::debug)
    {
        Pout<< FUNCTION_NAME << " : name:" << name << endl;
    }
    mode_t m = mode(name, checkVariants, followLink);

    if (S_ISREG(m))
    {
        return fileType::file;
    }
    else if (S_ISLNK(m))
    {
        return fileType::link;
    }
    else if (S_ISDIR(m))
    {
        return fileType::directory;
    }
    else
    {
        return fileType::undefined;
    }
}


bool Foam::exists
(
    const fileName& name,
    const bool checkVariants,
    const bool followLink
)
{
    if (POSIX::debug)
    {
        Pout<< FUNCTION_NAME << " : name:" << name << " checkVariants:"
            << bool(checkVariants) << " followLink:" << followLink << endl;
        if ((POSIX::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }
    return mode(name, checkVariants, followLink);
}


bool Foam::isDir(const fileName& name, const bool followLink)
{
    if (POSIX::debug)
    {
        Pout<< FUNCTION_NAME << " : name:" << name << " followLink:"
            << followLink << endl;
        if ((POSIX::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }
    return S_ISDIR(mode(name, false, followLink));
}


bool Foam::isFile
(
    const fileName& name,
    const bool checkVariants,
    const bool followLink
)
{
    if (POSIX::debug)
    {
        Pout<< FUNCTION_NAME << " : name:" << name << " checkVariants:"
            << bool(checkVariants) << " followLink:" << followLink << endl;
        if ((POSIX::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    return S_ISREG(mode(name, checkVariants, followLink));
}


off_t Foam::fileSize
(
    const fileName& name,
    const bool checkVariants,
    const bool followLink
)
{
    if (POSIX::debug)
    {
        Pout<< FUNCTION_NAME << " : name:" << name << " checkVariants:"
            << bool(checkVariants) << " followLink:" << followLink << endl;
        if ((POSIX::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }
    fileStat fileStatus(name, checkVariants, followLink);
    if (fileStatus.isValid())
    {
        return fileStatus.status().st_size;
    }
    else
    {
        return -1;
    }
}


time_t Foam::lastModified
(
    const fileName& name,
    const bool checkVariants,
    const bool followLink
)
{
    if (POSIX::debug)
    {
        Pout<< FUNCTION_NAME << " : name:" << name << " checkVariants:"
            << bool(checkVariants) << " followLink:" << followLink << endl;
        if ((POSIX::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }
    fileStat fileStatus(name, checkVariants, followLink);
    if (fileStatus.isValid())
    {
        return fileStatus.status().st_mtime;
    }
    else
    {
        return 0;
    }
}


double Foam::highResLastModified
(
    const fileName& name,
    const bool checkVariants,
    const bool followLink
)
{
    if (POSIX::debug)
    {
        Pout<< FUNCTION_NAME << " : name:" << name << " checkVariants:"
            << bool(checkVariants) << " followLink:" << followLink << endl;
        if ((POSIX::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }
    fileStat fileStatus(name, checkVariants, followLink);
    if (fileStatus.isValid())
    {
        return
            fileStatus.status().st_mtime
          + 1e-9*fileStatus.status().st_atim.tv_nsec;
    }
    else
    {
        return 0;
    }
}


Foam::fileNameList Foam::readDir
(
    const fileName& directory,
    const fileType type,
    const bool filterVariants,
    const bool followLink
)
{
    if (POSIX::debug)
    {
        Pout<< FUNCTION_NAME << " : reading directory " << directory << endl;
        if ((POSIX::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    // Create empty set of file names
    HashSet<fileName> dirEntries;

    // Pointers to the directory entries
    DIR *source;
    struct dirent *list;

    // Attempt to open directory and set the structure pointer
    if ((source = ::opendir(directory.c_str())) == nullptr)
    {
        if (POSIX::debug)
        {
            InfoInFunction
                << "cannot open directory " << directory << endl;
        }
    }
    else
    {
        // Read and parse all the entries in the directory
        while ((list = ::readdir(source)) != nullptr)
        {
            fileName fName(list->d_name);

            // Ignore files beginning with ., i.e. '.', '..' and '.*'
            if (fName.size() && fName[0] != '.')
            {
                word fExt = fName.ext();

                if
                (
                    (type == fileType::directory)
                 ||
                    (
                        type == fileType::file
                     && fName[fName.size()-1] != '~'
                     && fExt != "bak"
                     && fExt != "BAK"
                     && fExt != "old"
                     && fExt != "save"
                    )
                )
                {
                    if ((directory/fName).type(false, followLink) == type)
                    {
                        bool filtered = false;

                        if (filterVariants)
                        {
                            for (label i = 0; i < fileStat::nVariants_; ++ i)
                            {
                                if (fExt == fileStat::variantExts_[i])
                                {
                                    dirEntries.insert(fName.lessExt());
                                    filtered = true;
                                    break;
                                }
                            }
                        }

                        if (!filtered)
                        {
                            dirEntries.insert(fName);
                        }
                    }
                }
            }
        }

        ::closedir(source);
    }

    return dirEntries.toc();
}


bool Foam::cp(const fileName& src, const fileName& dest, const bool followLink)
{
    if (POSIX::debug)
    {
        Pout<< FUNCTION_NAME << " : src:" << src << " dest:" << dest << endl;
        if ((POSIX::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }
    // Make sure source exists.
    if (!exists(src))
    {
        return false;
    }

    const fileType srcType = src.type(false, followLink);

    fileName destFile(dest);

    // Check type of source file.
    if (srcType == fileType::file)
    {
        // If dest is a directory, create the destination file name.
        if (destFile.type() == fileType::directory)
        {
            destFile = destFile/src.name();
        }

        // Make sure the destination directory exists.
        if (!isDir(destFile.path()) && !mkDir(destFile.path()))
        {
            return false;
        }

        // Open and check streams.
        std::ifstream srcStream(src.c_str());
        if (!srcStream)
        {
            return false;
        }

        std::ofstream destStream(destFile.c_str());
        if (!destStream)
        {
            return false;
        }

        // Copy character data.
        destStream << srcStream.rdbuf();

        // Final check.
        if (!srcStream.eof() || !destStream)
        {
            return false;
        }
    }
    else if (srcType == fileType::link)
    {
        // If dest is a directory, create the destination file name.
        if (destFile.type() == fileType::directory)
        {
            destFile = destFile/src.name();
        }

        // Make sure the destination directory exists.
        if (!isDir(destFile.path()) && !mkDir(destFile.path()))
        {
            return false;
        }

        ln(src, destFile);
    }
    else if (srcType == fileType::directory)
    {
        // If dest is a directory, create the destination file name.
        if (destFile.type() == fileType::directory)
        {
            destFile = destFile/src.component(src.components().size() -1);
        }

        // Make sure the destination directory exists.
        if (!isDir(destFile) && !mkDir(destFile))
        {
            return false;
        }

        char* realSrcPath = realpath(src.c_str(), nullptr);
        char* realDestPath = realpath(destFile.c_str(), nullptr);
        const bool samePath = strcmp(realSrcPath, realDestPath) == 0;

        if (POSIX::debug && samePath)
        {
            InfoInFunction
                << "Attempt to copy " << realSrcPath << " to itself" << endl;
        }

        if (realSrcPath)
        {
            free(realSrcPath);
        }

        if (realDestPath)
        {
            free(realDestPath);
        }

        // Do not copy over self when src is actually a link to dest
        if (samePath)
        {
            return false;
        }

        // Copy files
        fileNameList contents = readDir(src, fileType::file, false, followLink);
        forAll(contents, i)
        {
            if (POSIX::debug)
            {
                InfoInFunction
                    << "Copying : " << src/contents[i]
                    << " to " << destFile/contents[i] << endl;
            }

            // File to file.
            cp(src/contents[i], destFile/contents[i], followLink);
        }

        // Copy sub directories.
        fileNameList subdirs = readDir
        (
            src,
            fileType::directory,
            false,
            followLink
        );

        forAll(subdirs, i)
        {
            if (POSIX::debug)
            {
                InfoInFunction
                    << "Copying : " << src/subdirs[i]
                    << " to " << destFile << endl;
            }

            // Dir to Dir.
            cp(src/subdirs[i], destFile, followLink);
        }
    }

    return true;
}


bool Foam::ln(const fileName& src, const fileName& dst)
{
    if (POSIX::debug)
    {
        Pout<< FUNCTION_NAME
            << " : Create softlink from : " << src << " to " << dst << endl;
        if ((POSIX::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    if (exists(dst))
    {
        WarningInFunction
            << "destination " << dst << " already exists. Not linking."
            << endl;
        return false;
    }

    if (src.isAbsolute() && !exists(src))
    {
        WarningInFunction
            << "source " << src << " does not exist." << endl;
        return false;
    }

    if (::symlink(src.c_str(), dst.c_str()) == 0)
    {
        return true;
    }
    else
    {
        WarningInFunction
            << "symlink from " << src << " to " << dst << " failed." << endl;
        return false;
    }
}


bool Foam::mv(const fileName& src, const fileName& dst, const bool followLink)
{
    if (POSIX::debug)
    {
        Pout<< FUNCTION_NAME << " : Move : " << src << " to " << dst << endl;
        if ((POSIX::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    if
    (
        dst.type() == fileType::directory
     && src.type(false, followLink) != fileType::directory
    )
    {
        const fileName dstName(dst/src.name());

        return ::rename(src.c_str(), dstName.c_str()) == 0;
    }
    else
    {
        return ::rename(src.c_str(), dst.c_str()) == 0;
    }
}


bool Foam::mvBak(const fileName& src, const std::string& ext)
{
    if (POSIX::debug)
    {
        Pout<< FUNCTION_NAME
            << " : moving : " << src << " to extension " << ext << endl;
        if ((POSIX::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    if (exists(src, false, false))
    {
        const int maxIndex = 99;
        char index[3];

        for (int n = 0; n <= maxIndex; n++)
        {
            fileName dstName(src + "." + ext);
            if (n)
            {
                sprintf(index, "%02d", n);
                dstName += index;
            }

            // Avoid overwriting existing files, except for the last
            // possible index where we have no choice
            if (!exists(dstName, false, false) || n == maxIndex)
            {
                return ::rename(src.c_str(), dstName.c_str()) == 0;
            }

        }
    }

    // Fall-through: nothing to do
    return false;
}


bool Foam::rm(const fileName& file)
{
    if (POSIX::debug)
    {
        Pout<< FUNCTION_NAME << " : Removing : " << file << endl;
        if ((POSIX::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    // Try returning plain file name; if not there, try variants
    if (remove(file.c_str()) == 0)
    {
        return true;
    }

    for (label i = 0; i < fileStat::nVariants_; ++ i)
    {
        const fileName fileVar = file + "." + fileStat::variantExts_[i];
        if (::remove(string(fileVar).c_str()) == 0)
        {
            return true;
        }
    }

    return false;
}


bool Foam::rmDir(const fileName& directory)
{
    if (POSIX::debug)
    {
        Pout<< FUNCTION_NAME << " : removing directory " << directory << endl;
        if ((POSIX::debug & 2) && !Pstream::master())
        {
            error::printStack(Pout);
        }
    }

    // Pointers to the directory entries
    DIR *source;
    struct dirent *list;

    // Attempt to open directory and set the structure pointer
    if ((source = ::opendir(directory.c_str())) == nullptr)
    {
        WarningInFunction
            << "cannot open directory " << directory << endl;

        return false;
    }
    else
    {
        // Read and parse all the entries in the directory
        while ((list = ::readdir(source)) != nullptr)
        {
            fileName fName(list->d_name);

            if (fName != "." && fName != "..")
            {
                fileName path = directory/fName;

                if (path.type(false, false) == fileType::directory)
                {
                    if (!rmDir(path))
                    {
                        WarningInFunction
                            << "failed to remove directory " << fName
                            << " while removing directory " << directory
                            << endl;

                        ::closedir(source);

                        return false;
                    }
                }
                else
                {
                    if (!rm(path))
                    {
                        WarningInFunction
                            << "failed to remove file " << fName
                            << " while removing directory " << directory
                            << endl;

                        ::closedir(source);

                        return false;
                    }
                }
            }

        }

        if (!rm(directory))
        {
            WarningInFunction
                << "failed to remove directory " << directory << endl;

            ::closedir(source);

            return false;
        }

        ::closedir(source);

        return true;
    }
}


unsigned int Foam::sleep(const unsigned int s)
{
    return ::sleep(s);
}


void Foam::fdClose(const int fd)
{
    if (close(fd) != 0)
    {
        FatalErrorInFunction
            << "close error on " << fd << endl
            << abort(FatalError);
    }
}


bool Foam::ping
(
    const string& destName,
    const label destPort,
    const label timeOut
)
{
    struct hostent *hostPtr;
    volatile int sockfd;
    struct sockaddr_in destAddr;      // Will hold the destination addr
    u_int addr;

    if ((hostPtr = ::gethostbyname(destName.c_str())) == nullptr)
    {
        FatalErrorInFunction
            << "gethostbyname error " << h_errno << " for host " << destName
            << abort(FatalError);
    }

    // Get first of the SLL of addresses
    addr = (reinterpret_cast<struct in_addr*>(*(hostPtr->h_addr_list)))->s_addr;

    // Allocate socket
    sockfd = ::socket(AF_INET, SOCK_STREAM, 0);
    if (sockfd < 0)
    {
        FatalErrorInFunction
            << "socket error"
            << abort(FatalError);
    }

    // Fill sockaddr_in structure with dest address and port
    memset(reinterpret_cast<char *>(&destAddr), '\0', sizeof(destAddr));
    destAddr.sin_family = AF_INET;
    destAddr.sin_port = htons(ushort(destPort));
    destAddr.sin_addr.s_addr = addr;


    timer myTimer(timeOut);

    if (timedOut(myTimer))
    {
        // Setjmp from timer jumps back to here
        fdClose(sockfd);
        return false;
    }

    if
    (
        ::connect
        (
            sockfd,
            reinterpret_cast<struct sockaddr*>(&destAddr),
            sizeof(struct sockaddr)
        ) != 0
    )
    {
        // Connection refused. Check if network was actually used or not.

        int connectErr = errno;

        fdClose(sockfd);

        if (connectErr == ECONNREFUSED)
        {
            return true;
        }

        return false;
    }

    fdClose(sockfd);

    return true;
}


bool Foam::ping(const string& hostname, const label timeOut)
{
    return ping(hostname, 222, timeOut) || ping(hostname, 22, timeOut);
}


int Foam::system(const std::string& command)
{
    return ::system(command.c_str());
}


void* Foam::dlOpen(const fileName& lib, const bool check)
{
    if (POSIX::debug)
    {
        std::cout<< "dlOpen(const fileName&)"
            << " : dlopen of " << lib << std::endl;
    }
    void* handle = ::dlopen(lib.c_str(), RTLD_LAZY|RTLD_GLOBAL);

    if (!handle && check)
    {
        WarningInFunction
            << "dlopen error : " << ::dlerror()
            << endl;
    }

    if (POSIX::debug)
    {
        std::cout
            << "dlOpen(const fileName&)"
            << " : dlopen of " << lib
            << " handle " << handle << std::endl;
    }

    return handle;
}


bool Foam::dlClose(void* handle)
{
    if (POSIX::debug)
    {
        std::cout
            << "dlClose(void*)"
            << " : dlclose of handle " << handle << std::endl;
    }
    return ::dlclose(handle) == 0;
}


void* Foam::dlSym(void* handle, const std::string& symbol)
{
    if (POSIX::debug)
    {
        std::cout
            << "dlSym(void*, const std::string&)"
            << " : dlsym of " << symbol << std::endl;
    }

    // Clear any old errors - see manpage dlopen
    (void) ::dlerror();

    // Get address of symbol
    void* fun = ::dlsym(handle, symbol.c_str());

    // Find error (if any)
    char *error = ::dlerror();

    if (error)
    {
        WarningInFunction
            << "Cannot lookup symbol " << symbol << " : " << error
            << endl;
    }

    return fun;
}


bool Foam::dlSymFound(void* handle, const std::string& symbol)
{
    if (handle && !symbol.empty())
    {
        if (POSIX::debug)
        {
            std::cout
                << "dlSymFound(void*, const std::string&)"
                << " : dlsym of " << symbol << std::endl;
        }

        // Clear any old errors - see manpage dlopen
        (void) ::dlerror();

        // Get address of symbol
        (void) ::dlsym(handle, symbol.c_str());

        // Symbol can be found if there was no error
        return !::dlerror();
    }
    else
    {
        return false;
    }
}


static int collectLibsCallback
(
    struct dl_phdr_info *info,
    size_t size,
    void *data
)
{
    Foam::DynamicList<Foam::fileName>* ptr =
        reinterpret_cast<Foam::DynamicList<Foam::fileName>*>(data);
    ptr->append(info->dlpi_name);
    return 0;
}


Foam::fileNameList Foam::dlLoaded()
{
    DynamicList<fileName> libs;
    dl_iterate_phdr(collectLibsCallback, &libs);
    if (POSIX::debug)
    {
        std::cout
            << "dlLoaded()"
            << " : determined loaded libraries :" << libs.size() << std::endl;
    }

    return move(libs);
}


// ************************************************************************* //
