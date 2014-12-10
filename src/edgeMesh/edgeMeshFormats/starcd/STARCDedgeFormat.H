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
    Foam::fileFormats::STARCDedgeFormat

Description
    Read/write the lines from pro-STAR vrt/cel files.

Note
    Uses the extension \a .inp (input) to denote the format.

See Also
    Foam::meshReaders::STARCD

SourceFiles
    STARCDedgeFormat.C

\*---------------------------------------------------------------------------*/

#ifndef STARCDedgeFormat_H
#define STARCDedgeFormat_H

#include "edgeMesh.H"
#include "STARCDCore.H"

#include "IFstream.H"
#include "Ostream.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                      Class STARCDedgeFormat Declaration
\*---------------------------------------------------------------------------*/

class STARCDedgeFormat
:
    public edgeMesh,
    public STARCDCore
{
    // Private Data

        //- STAR-CD identifier for line shapes (1d elements)
        static const int starcdLineShape_ = 2;

        //- STAR-CD identifier for line type
        static const int starcdLineType_  = 5;


    // Private Member Functions

        static inline void writeLines
        (
            Ostream&,
            const edgeList&
        );

        //- Disallow default bitwise copy construct
        STARCDedgeFormat(const STARCDedgeFormat&);

        //- Disallow default bitwise assignment
        void operator=(const STARCDedgeFormat&);


protected:

    // Protected Member Functions

    static void writeCase
    (
        Ostream&,
        const pointField&,
        const label nEdges
    );

public:

    // Constructors

        //- Construct from file name
        STARCDedgeFormat(const fileName&);


    // Selectors

        //- Read file and return edgeMesh
        static autoPtr<edgeMesh> New(const fileName& name)
        {
            return autoPtr<edgeMesh>
            (
                new STARCDedgeFormat(name)
            );
        }


    //- Destructor
    virtual ~STARCDedgeFormat()
    {}


    // Member Functions

        //- Write edge mesh
        static void write(const fileName&, const edgeMesh&);

        //- Read from file
        virtual bool read(const fileName&);

        //- Write object
        virtual void write(const fileName& name) const
        {
            write(name, *this);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
