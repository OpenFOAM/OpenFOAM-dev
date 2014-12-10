/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::ensightFile

Description
    Ensight output with specialized write() for strings, integers and floats.
    Correctly handles binary write as well.

\*---------------------------------------------------------------------------*/

#ifndef ensightFile_H
#define ensightFile_H

#include "OFstream.H"
#include "IOstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class ensightFile Declaration
\*---------------------------------------------------------------------------*/

class ensightFile
:
    public OFstream
{
    // Private data

        //- allow undef in results
        static bool allowUndef_;

        //- value to represent undef in results
        static scalar undefValue_;

        //- The '*' mask appropriate for subDir
        static string mask_;

        //- The printf format for zero-padded subdirectory numbers
        static string dirFmt_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const ensightFile&);

        //- Disallow default copy constructor
        ensightFile(const ensightFile&);


public:

    // Constructors

        //- Construct from pathname
        ensightFile
        (
            const fileName& pathname,
            IOstream::streamFormat format=IOstream::BINARY
        );


    //- Destructor
    ~ensightFile();


    // Access

        //- Return setting for whether 'undef' values are allowed in results
        static bool allowUndef();

        //- The '*' mask appropriate for subDir
        static string mask();

        //- Consistent zero-padded numbers for subdirectories
        static string subDir(const label);

        //- Set width of subDir and mask. Default width is 8 digits.
        //  Max width is 31 digits.
        static void subDirWidth(const label);

        //- Return current width of subDir and mask.
        static label subDirWidth();


    // Edit

        static bool allowUndef(bool);

        //- Assign the value to represent undef in the results
        //  Returns the previous value
        //  NB: do not use values larger than floatScalarVGREAT
        static scalar undefValue(const scalar);


    // Output

        //- binary write
        virtual Ostream& write(const char* buf, std::streamsize count);

        //- write element keyword with trailing newline, optionally with undef
        virtual Ostream& writeKeyword(const string& key);

        //- write "C Binary" for binary files (eg, geometry/measured)
        Ostream& writeBinaryHeader();

        //- write undef value
        Ostream& writeUndef();

        //- write string as "%80s" or as binary
        Ostream& write(const string& value);

        //- write integer as "%10d" or as binary
        Ostream& write(const label value);

        //- write integer with specified width or as binary
        Ostream& write(const label value, const label fieldWidth);

        //- write float as "%12.5e" or as binary
        Ostream& write(const scalar value);

        //- Add carriage return to ascii stream
        void newline();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
