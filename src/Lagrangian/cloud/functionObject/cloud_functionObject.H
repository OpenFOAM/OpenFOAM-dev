/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::functionObjects::cloud

Description
    Function object that solves for the evolution of a cloud. Only provides
    one-way coupling with a finite-volume carrier phase. The carrier phase
    affects the cloud, but the cloud does not affect the carrier phase.

SourceFiles
    cloud_functionObject.C
    cloud_functionObjectTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef cloud_functionObject_H
#define cloud_functionObject_H

#include "regionFunctionObject.H"
#include "cloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                            Class cloud Declaration
\*---------------------------------------------------------------------------*/

class cloud
:
    public regionFunctionObject
{
    // Private Data

        //- The cloud
        autoPtr<Foam::cloud> cloudPtr_;


protected:

    // Protected Classes

        //- Empty template class used to template the constructor on the type
        //  of the cloud
        template<class Type>
        class Cloud
        {};


    // Protected Constructors

        //- Construct from components
        template<class Type>
        cloud
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict,
            const Cloud<Type>&
        );


public:

    // Public Static Data

        //- Run-time type information
        TypeName("cloud");


    //- Destructor
    virtual ~cloud();


    // Member Functions

        //- Return the list of fields required
        virtual wordList fields() const;

        //- Return false so that this function does not execute at the start
        virtual bool executeAtStart() const;

        //- Solve the cloud
        virtual bool execute();

        //- Dummy write
        virtual bool write();


        // Mesh changes

            //- Update for mesh motion
            virtual void movePoints(const polyMesh&);

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&);
};


/*---------------------------------------------------------------------------*\
                            Class Cloud Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Cloud
:
    public cloud
{
public:

    // Public Static Data

        //- Run-time type information
        TypeName("Cloud");


    // Constructors

        //- Construct from components
        Cloud
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "cloud_functionObjectTemplates.C"
#endif


#define makeCloudFunctionObject(Type)                                       \
                                                                            \
    typedef Cloud<clouds::Type> Type##CloudFunctionObject;                  \
                                                                            \
    defineTemplateTypeNameAndDebugWithName                                  \
    (                                                                       \
        Type##CloudFunctionObject,                                          \
        (clouds::Type::typeName + "Cloud").c_str(),                         \
        0                                                                   \
    );                                                                      \
                                                                            \
    addToRunTimeSelectionTable                                              \
    (                                                                       \
        functionObject,                                                     \
        Type##CloudFunctionObject,                                          \
        dictionary                                                          \
    );


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
