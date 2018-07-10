/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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
    Foam::functionObjects::sixDoFRigidBodyState

Description
    Writes the 6-DoF motion state.

    Example of function object specification:
    \verbatim
    sixDoFRigidBodyState
    {
        type           sixDoFRigidBodyState;
        libs           ("libsixDoFRigidBodyState.so");
        angleFormat    degrees;
    }
    \endverbatim

Usage
    \table
        Property     | Description                  | Required | Default value
        type         | type name: sixDoFRigidBodyState    | yes |
        angleFormat  | degrees or radians           | no       | radian
    \endtable

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::logFiles

SourceFiles
    sixDoFRigidBodyState.C

\*---------------------------------------------------------------------------*/

#ifndef sixDoFRigidBodyState_H
#define sixDoFRigidBodyState_H

#include "fvMeshFunctionObject.H"
#include "logFiles.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class sixDoFRigidBodyState Declaration
\*---------------------------------------------------------------------------*/

class sixDoFRigidBodyState
:
    public fvMeshFunctionObject,
    public logFiles
{
    // Private data

        word angleFormat_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        sixDoFRigidBodyState(const sixDoFRigidBodyState&);

        //- Disallow default bitwise assignment
        void operator=(const sixDoFRigidBodyState&);


protected:

    // Protected Member Functions

        //- overloaded writeFileHeader from writeFile
        virtual void writeFileHeader(const label i = 0);


public:

    //- Runtime type information
    TypeName("sixDoFRigidBodyState");


    // Constructors

        //- Construct from Time and dictionary
        sixDoFRigidBodyState
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~sixDoFRigidBodyState();


    // Member Functions

        //- Read the sixDoFRigidBodyState data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the sixDoFRigidBodyState
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
