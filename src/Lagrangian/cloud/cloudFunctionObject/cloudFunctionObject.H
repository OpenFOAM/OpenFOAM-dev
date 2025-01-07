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
    Foam::functionObjects::cloudFunctionObject

Description
    Base class for function objects that refer to a cloud. Provides hooks into
    the cloud solution process. Not a function object in itself. Must be mixed
    in.

SourceFiles
    cloudFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_cloudFunctionObject_H
#define functionObjects_cloudFunctionObject_H

#include "LagrangianMeshFunctionObject.H"
#include "fvMeshFunctionObject.H"
#include "LagrangianFieldsFwd.H"
#include "LagrangianSubFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class cloud;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class cloudFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class cloudFunctionObject
{
private:

    // Private Data

        //- Reference to the function object
        const functionObject& function_;

        //- Reference to the cloud
        const Foam::cloud& cloud_;


protected:

    // Protected Member Functions

        //- Determine whether the cloud is of the given type
        template<class Cloud>
        bool isCloud() const;

        //- Get a reference to the cloud
        template<class Cloud = Foam::cloud>
        const Cloud& cloud() const;


public:

    //- Runtime type information
    TypeName("cloudFunctionObject");


    // Constructors

        //- Construct from a Lagrangian mesh function object
        cloudFunctionObject(const LagrangianMeshFunctionObject&);

        //- Construct from an fvMesh function object and a dictionary
        cloudFunctionObject
        (
            const fvMeshFunctionObject&,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        cloudFunctionObject
        (
            const cloudFunctionObject&
        ) = delete;


    //- Destructor
    virtual ~cloudFunctionObject();


    // Member Functions

        //- Return the cloud
        const Foam::cloud& cloud() const;

        //- Hook before solution steps
        virtual void preSolve();

        //- Solve equations and/or update continually changing properties
        virtual void calculate
        (
            const LagrangianSubScalarField& deltaT,
            const bool final
        );

        //- Hook before all face crossings
        virtual void preCrossFaces
        (
            const LagrangianScalarInternalDynamicField& fraction
        );

        //- Hook before face crossings of a specific sub-mesh
        virtual void preCrossFaces
        (
            const LagrangianSubScalarSubField& fraction
        );

        //- Hook following face crossings of a specific sub-mesh
        virtual void postCrossFaces
        (
            const LagrangianSubScalarSubField& fraction
        );

        //- Hook after all face crossings
        virtual void postCrossFaces
        (
            const LagrangianScalarInternalDynamicField& fraction
        );

        //- Hook following solution steps
        virtual void postSolve();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const cloudFunctionObject&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "cloudFunctionObjectTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
