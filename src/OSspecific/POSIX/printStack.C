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

#include "error.H"
#include "IStringStream.H"
#include "OStringStream.H"
#include "OSspecific.H"
#include "IFstream.H"
#include "ReadHex.H"

#include <cxxabi.h>
#include <execinfo.h>
#include <dlfcn.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

string pOpen(const string &cmd, label line=0)
{
    const int MAX = 1000;

    FILE *cmdPipe = popen(cmd.c_str(), "r");

    if (cmdPipe)
    {
        // Read line number of lines
        for (label cnt = 0; cnt <= line; cnt++)
        {
            char buffer[MAX];
            char* s = fgets(buffer, MAX-1, cmdPipe);

            if (s == NULL)
            {
                return "";
            }

            if (cnt == line)
            {
                string str(buffer);
                return str.substr(0, str.size()-1);
            }
        }
        pclose(cmdPipe);
    }

    return "";
}


// use popen to call addr2line (using bfd.h directly would have
// meant relinking everything)

void printSourceFileAndLine
(
    Ostream& os,
    const HashTable<label, fileName>& addressMap,
    const fileName& filename,
    const word& address
)
{
    word myAddress = address;

    if (filename.ext() == "so")
    {
        // Convert offset into .so into offset into executable.

        void *addr;
        sscanf(myAddress.c_str(), "%p",&addr);

        Dl_info info;

        dladdr(addr, &info);

        unsigned long offset = ulong(info.dli_fbase);

        IStringStream addressStr(address.substr(2));
        long addressValue = ReadHex<long>(addressStr);
        long relativeAddress = addressValue-offset;

        // Reconstruct hex word from address
        OStringStream nStream;
        nStream << "0x" << hex << relativeAddress;
        myAddress = nStream.str();
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


void getSymbolForRaw
(
    Ostream& os,
    const string& raw,
    const fileName& filename,
    const word& address
)
{
    if (filename.size() && filename[0] == '/')
    {
        string fcnt = pOpen
        (
            "addr2line -f --demangle=auto --exe "
          + filename
          + " "
          + address
        );

        if (fcnt != "")
        {
            os  << fcnt.c_str();
            return;
        }
    }
    os  << "Uninterpreted: " << raw.c_str();
}


void error::safePrintStack(std::ostream& os)
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


void error::printStack(Ostream& os)
{
    // Reads the starting addresses for the dynamically linked libraries
    // from the /proc/pid/maps-file
    // I'm afraid this works only for Linux 2.6-Kernels (may work on 2.4)
    // Note2: the filenames in here will have softlinks resolved so will
    // go wrong when having e.g. OpenFOAM installed under a softlink.

    HashTable<label, fileName> addressMap;
    {
        IFstream is("/proc/" + name(pid()) + "/maps");

        while (is.good())
        {
            string line;
            is.getLine(line);

            string::size_type space = line.rfind(' ') + 1;
            fileName libPath = line.substr(space, line.size()-space);

            if (libPath.size() && libPath[0] == '/')
            {
                string offsetString(line.substr(0, line.find('-')));
                IStringStream offsetStr(offsetString);
                addressMap.insert(libPath, ReadHex<label>(offsetStr));
            }
        }
    }

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

        os  << '#' << label(i) << "  ";
        //os  << "Raw   : " << msg << "\n\t";
        {
            string::size_type lPos = msg.find('[');
            string::size_type rPos = msg.find(']');

            if (lPos != string::npos && rPos != string::npos && lPos < rPos)
            {
                address = msg.substr(lPos+1, rPos-lPos-1);
                msg = msg.substr(0, lPos);
            }

            string::size_type bracketPos = msg.find('(');
            string::size_type spacePos = msg.find(' ');
            if (bracketPos != string::npos || spacePos != string::npos)
            {
                programFile = msg.substr(0, min(spacePos, bracketPos));

                // not an absolute path
                if (programFile[0] != '/')
                {
                    string tmp = pOpen("which " + programFile);
                    if (tmp[0] == '/' || tmp[0] == '~')
                    {
                        programFile = tmp;
                    }
                }
            }
        }

        string::size_type bracketPos = msg.find('(');

        if (bracketPos != string::npos)
        {
            string::size_type start = bracketPos+1;

            string::size_type plusPos = msg.find('+', start);

            if (plusPos != string::npos)
            {
                string cName(msg.substr(start, plusPos-start));

                int status;
                char* cplusNamePtr = abi::__cxa_demangle
                (
                    cName.c_str(),
                    NULL,                   // have it malloc itself
                    0,
                    &status
                );

                if (status == 0 && cplusNamePtr)
                {
                    os  << cplusNamePtr;
                    free(cplusNamePtr);
                }
                else
                {
                    os  << cName.c_str();
                }
            }
            else
            {
                string::size_type endBracketPos = msg.find(')', start);

                if (endBracketPos != string::npos)
                {
                    string fullName(msg.substr(start, endBracketPos-start));

                    os  << fullName.c_str() << nl;
                }
                else
                {
                    // Print raw message
                    getSymbolForRaw(os, msg, programFile, address);
                }
            }
        }
        else
        {
            // Print raw message
            getSymbolForRaw(os, msg, programFile, address);
        }

        printSourceFileAndLine(os, addressMap, programFile, address);

        os  << nl;
    }

    free(strings);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
