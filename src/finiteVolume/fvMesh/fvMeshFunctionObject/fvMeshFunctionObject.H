/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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
    Foam::functionObjects::fvMeshFunctionObject

Description
    Specialization of Foam::functionObject for an Foam::fvMesh, providing a
    reference to the Foam::fvMesh.

    If the selected region is not an Foam::fvMesh a Foam::FatalError will be
    generated.

SeeAlso
    Foam::regionFunctionObject
    Foam::functionObject

SourceFiles
    fvMeshFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_fvMeshFunctionObject_H
#define functionObjects_fvMeshFunctionObject_H

#include "regionFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                     Class fvMeshFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class fvMeshFunctionObject
:
    public regionFunctionObject
{

protected:

    // Protected member data

        //- Reference to the fvMesh
        const fvMesh& mesh_;


    // Protected member functions

        //- Find field in the objectRegistry
        template<class FieldType>
        bool foundField(const word& fieldName) const;

        //- Lookup field from the objectRegistry
        template<class FieldType>
        const FieldType& lookupField(const word& fieldName) const;

        //- Store the given field in the objectRegistry under the given name
        template<class FieldType>
        bool store
        (
            word& fieldName,
            const tmp<FieldType>& tfield,
            bool cacheable = false
        );

        //- Write field if present in objectRegistry
        bool writeField(const word& fieldName);


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        fvMeshFunctionObject(const fvMeshFunctionObject&);

        //- Disallow default bitwise assignment
        void operator=(const fvMeshFunctionObject&);


public:

    //- Runtime type information
    TypeName("fvMeshFunctionObject");


    // Constructors

        //- Construct from Time and dictionary
        fvMeshFunctionObject
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~fvMeshFunctionObject();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvMeshFunctionObjectTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
