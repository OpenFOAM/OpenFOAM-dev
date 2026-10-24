/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::preLoad

Description
    preLoad provides a mechanism to create and maintain lists of dynamically
    generated and cached objects which need to be pre-loaded on restart,
    before any mesh changes occur.

    The preLoad class is a form of mesh object created by the preload::New
    method and stored on the fvMesh.  The names of fields which need to be
    pre-loaded on restart are added directly to the preLoad::fields
    HashSet<word> and this list is written out automatically into the
    <time>/uniform/preLoad dictionary.  On restart the constructor of the
    preLoad object reads the <time>/uniform/preLoad dictionary and automatically
    constructs and registers the list of fields specified in the
    \c fields entry.

SourceFiles
    preLoad.C

\*---------------------------------------------------------------------------*/

#ifndef preLoad_H
#define preLoad_H

#include "timeIOdictionary.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                          Class preLoad Declaration
\*---------------------------------------------------------------------------*/

class preLoad
:
    public timeIOdictionary
{
    // Protected Constructors

        //- Construct from fvMesh
        explicit preLoad(const fvMesh& mesh);


public:

    //- Runtime type information
    TypeName("preLoad");

    //- Set of field names to pre-load
    HashSet<word> fields;


    // Constructors

        //- Disallow default bitwise copy construction
        preLoad(const preLoad&) = delete;

        //- Construct and return the preLoad object
        static preLoad& New(const fvMesh& mesh);


    //- Destructor
    virtual ~preLoad()
    {}


    // Member Functions

        //- Write using given format, version and compression
        //  but only if fields is not empty
        virtual bool writeObject
        (
            IOstream::streamFormat,
            IOstream::versionNumber,
            IOstream::compressionType,
            const bool updateInstance
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const preLoad&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
