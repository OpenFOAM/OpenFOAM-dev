/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::fileFormats::extendedEdgeMeshFormat

Description
    Provide a means of reading/writing the single-file OpenFOAM
    extendedEdgeMesh format

SourceFiles
    extendedEdgeMeshFormat.C

\*---------------------------------------------------------------------------*/

#ifndef extendedEdgeMeshFormat_H
#define extendedEdgeMeshFormat_H

#include "extendedEdgeMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                        Class extendedEdgeMeshFormat Declaration
\*---------------------------------------------------------------------------*/

class extendedEdgeMeshFormat
:
    public extendedEdgeMesh
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        extendedEdgeMeshFormat(const extendedEdgeMeshFormat&);

        //- Disallow default bitwise assignment
        void operator=(const extendedEdgeMeshFormat&);


public:

    // Constructors

        //- Construct from file name
        extendedEdgeMeshFormat(const fileName&);


    //- Destructor
    virtual ~extendedEdgeMeshFormat()
    {}


    // Member Functions

        //- Read from file
        virtual bool read(const fileName&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
