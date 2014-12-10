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
    Foam::fileFormats::edgeMeshFormatsCore

Description
    A collection of helper functions for reading/writing edge formats.

SourceFiles
    edgeMeshFormatsCore.C

\*---------------------------------------------------------------------------*/

#ifndef edgeMeshFormatsCore_H
#define edgeMeshFormatsCore_H

#include "Map.H"
#include "HashSet.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

class IFstream;
class Time;

namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                       Class edgeMeshFormatsCore Declaration
\*---------------------------------------------------------------------------*/

class edgeMeshFormatsCore
{
protected:

    // Protected Member Functions

        //- Read non-comment line
        static string getLineNoComment(IFstream&);

public:

    // Static Data

        //- The file extension corresponding to 'native' edge format
        //  Normally "eMesh" (edge-mesh)
        static word nativeExt;

    // Static Member Functions

        static bool checkSupport
        (
            const wordHashSet& available,
            const word& ext,
            const bool verbose,
            const word& functionName
        );

       // //- Return the local file name (within time directory)
       // //  NEEDS FIXING
       // static fileName localMeshFileName(const word& edgeName="");

       // //- Find instance with edgeName
       // //  NEEDS FIXING
       // static fileName findMeshInstance
       // (
       //     const Time&,
       //     const word& edgeName=""
       // );

       // //- Find mesh file with edgeName
       // //  NEEDS FIXING
       // static fileName findMeshFile(const Time&, const word& edgeName="");


    // Constructors

        //- Construct null
        edgeMeshFormatsCore();


    //- Destructor
    virtual ~edgeMeshFormatsCore();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
