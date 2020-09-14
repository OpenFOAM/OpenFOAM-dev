/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::OFstream

Description
    Output to file stream.

SourceFiles
    OFstream.C

\*---------------------------------------------------------------------------*/

#ifndef OFstream_H
#define OFstream_H

#include "OSstream.H"
#include "fileName.H"
#include "className.H"

#include <fstream>
using std::ofstream;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class OFstream;

/*---------------------------------------------------------------------------*\
                      Class OFstreamAllocator Declaration
\*---------------------------------------------------------------------------*/

//- A std::ostream with ability to handle compressed files
class OFstreamAllocator
{
    friend class OFstream;

    ostream* ofPtr_;

    // Constructors

        //- Construct from filePath
        OFstreamAllocator
        (
            const fileName& filePath,
            IOstream::compressionType compression=IOstream::UNCOMPRESSED,
            const bool append = false
        );

    //- Destructor
    ~OFstreamAllocator();

};


/*---------------------------------------------------------------------------*\
                          Class OFstream Declaration
\*---------------------------------------------------------------------------*/

class OFstream
:
    private OFstreamAllocator,
    public OSstream
{
    // Private Data

        fileName filePath_;


public:

    // Declare name of the class and its debug switch
    ClassName("OFstream");


    // Constructors

        //- Construct from filePath
        OFstream
        (
            const fileName& filePath,
            streamFormat format=ASCII,
            versionNumber version=currentVersion,
            compressionType compression=UNCOMPRESSED,
            const bool append = false
        );


    //- Destructor
    ~OFstream();


    // Member Functions

        // Access

            //- Return the name of the stream
            const fileName& name() const
            {
                return filePath_;
            }

            //- Return non-const access to the name of the stream
            fileName& name()
            {
                return filePath_;
            }


        // STL stream

            //- Access to underlying std::ostream
            virtual ostream& stdStream();

            //- Const access to underlying std::ostream
            virtual const ostream& stdStream() const;


        // Print

            //- Print description of IOstream to Ostream
            void print(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Global predefined null output stream "/dev/null"
extern OFstream Snull;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
