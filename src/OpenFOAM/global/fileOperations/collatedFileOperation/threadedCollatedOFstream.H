/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2024 OpenFOAM Foundation
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
    Foam::threadedCollatedOFstream

Description
    Master-only drop-in replacement for OFstream.

SourceFiles
    threadedCollatedOFstream.C

\*---------------------------------------------------------------------------*/

#ifndef threadedCollatedOFstream_H
#define threadedCollatedOFstream_H

#include "OStringStream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class OFstreamCollator;

/*---------------------------------------------------------------------------*\
                  Class threadedCollatedOFstream Declaration
\*---------------------------------------------------------------------------*/

class threadedCollatedOFstream
:
    public OStringStream
{
    // Private Data

        OFstreamCollator& writer_;

        const fileName filePath_;

        const IOstream::compressionType compression_;

        const bool useThread_;


public:

    // Constructors

        //- Construct and set stream status
        threadedCollatedOFstream
        (
            OFstreamCollator&,
            const fileName& filePath,
            const streamFormat format = ASCII,
            const versionNumber version = currentVersion,
            const compressionType compression = UNCOMPRESSED,
            const bool useThread = true
        );


    //- Destructor
    ~threadedCollatedOFstream();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
