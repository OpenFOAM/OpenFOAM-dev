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

#include "error.H"
#include "OStringStream.H"
#include "OSspecific.H"
#include "IFstream.H"

#include <inttypes.h>
#include <cxxabi.h>
#include <execinfo.h>
#include <dlfcn.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

string pOpen(const string &cmd, label line=0)
{
    string res = "\n";

    FILE *cmdPipe = popen(cmd.c_str(), "r");

    if (cmdPipe)
    {
        char *buf = nullptr;

        // Read line number of lines
        for (label cnt = 0; cnt <= line; cnt++)
        {
            size_t linecap = 0;
            ssize_t linelen;
            linelen = getline(&buf, &linecap, cmdPipe);

            if (linelen < 0)
            {
                break;
            }

            if (cnt == line)
            {
                res = string(buf);
                break;
            }
        }

        if (buf != nullptr)
        {
            free(buf);
        }

        pclose(cmdPipe);
    }

    return res.substr(0, res.size() - 1);
}


inline word addressToWord(const uintptr_t addr)
{
    OStringStream nStream;
    nStream << "0x" << hex << addr;
    return nStream.str();
}


void printSourceFileAndLine
(
    Ostream& os,
    const fileName& filename,
    Dl_info *info,
    void *addr
)
{
    uintptr_t address = uintptr_t(addr);
    word myAddress = addressToWord(address);

    if (filename.ext() == "so")
    {
        // Convert address into offset into dynamic library
        uintptr_t offset = uintptr_t(info->dli_fbase);
        intptr_t relativeAddress = address - offset;
        myAddress = addressToWord(relativeAddress);
    }

    if (filename[0] == '/')
    {
        string line = pOpen
        (
            "addr2line -f --demangle=auto --exe "
          + filename
          + " "
          + myAddress,
            1
        );

        if (line == "")
        {
            os  << " addr2line failed";
        }
        else if (line == "??:0")
        {
            os  << " in " << filename;
        }
        else
        {
            string cwdLine(line.replaceAll(cwd() + '/', ""));
            string homeLine(cwdLine.replaceAll(home(), '~'));

            os  << " at " << homeLine.c_str();
        }
    }
}


fileName absolutePath(const char* fn)
{
    fileName fname(fn);

    if (fname[0] != '/' && fname[0] != '~')
    {
        string tmp = pOpen("which " + fname);

        if (tmp[0] == '/' || tmp[0] == '~')
        {
            fname = tmp;
        }
    }

    return fname;
}


word demangleSymbol(const char* sn)
{
    word res;
    int st;
    char* cxx_sname = abi::__cxa_demangle
    (
        sn,
        nullptr,
        0,
        &st
    );

    if (st == 0 && cxx_sname)
    {
        res = word(cxx_sname);
        free(cxx_sname);
    }
    else
    {
        res = word(sn);
    }

    return res;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


void Foam::error::safePrintStack(std::ostream& os)
{
    // Get raw stack symbols
    void *array[100];
    size_t size = backtrace(array, 100);
    char **strings = backtrace_symbols(array, size);

    // See if they contain function between () e.g. "(__libc_start_main+0xd0)"
    // and see if cplus_demangle can make sense of part before +
    for (size_t i = 0; i < size; i++)
    {
        string msg(strings[i]);
        fileName programFile;
        word address;

        os  << '#' << label(i) << '\t' << msg << std::endl;
    }
}


void Foam::error::printStack(Ostream& os)
{
    // Get raw stack symbols
    const size_t CALLSTACK_SIZE = 128;

    void *callstack[CALLSTACK_SIZE];
    size_t size = backtrace(callstack, CALLSTACK_SIZE);

    Dl_info *info = new Dl_info;

    fileName fname = "???";
    word address;

    for(size_t i=0; i<size; i++)
    {
        int st = dladdr(callstack[i], info);

        os << '#' << label(i) << "  ";
        if (st != 0 && info->dli_fname != nullptr && info->dli_fname[0] != '\0')
        {
            fname = absolutePath(info->dli_fname);

            os <<
            (
                (info->dli_sname != nullptr)
              ? demangleSymbol(info->dli_sname)
              : "?"
            );
        }
        else
        {
            os << "?";
        }

        printSourceFileAndLine(os, fname, info, callstack[i]);
        os << nl;
    }

    delete info;
}


// ************************************************************************* //
