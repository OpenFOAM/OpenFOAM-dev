/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2020 OpenFOAM Foundation
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
    Foam::masterOFstream

Description
    Master-only drop-in replacement for OFstream.

SourceFiles
    masterOFstream.C

\*---------------------------------------------------------------------------*/

#ifndef masterOFstream_H
#define masterOFstream_H

#include "OStringStream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class masterOFstream Declaration
\*---------------------------------------------------------------------------*/

class masterOFstream
:
    public OStringStream
{
    // Private Data

        const fileName filePath_;

        const IOstream::compressionType compression_;

        const bool append_;

        //- Should file be written
        const bool write_;


    // Private Member Functions

        //- Open file with checking
        void checkWrite(const fileName& fName, const string& str);


public:

    // Constructors

        //- Construct and set stream status
        masterOFstream
        (
            const fileName& filePath,
            streamFormat format=ASCII,
            versionNumber version=currentVersion,
            compressionType compression=UNCOMPRESSED,
            const bool append = false,
            const bool write = true
        );


    //- Destructor
    ~masterOFstream();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
