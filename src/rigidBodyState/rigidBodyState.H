/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2020 OpenFOAM Foundation
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
    Foam::functionObjects::rigidBodyState

Description
    Writes the rigid body motion state.

    Example of function object specification:
    \verbatim
    rigidBodyState
    {
        type           rigidBodyState;
        libs           ("librigidBodyState.so");
        angleFormat    degrees;
    }
    \endverbatim

Usage
    \table
        Property     | Description                  | Required | Default value
        type         | type name: rigidBodyState    | yes |
        angleFormat  | degrees or radians           | no       | radian
    \endtable

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::logFiles

SourceFiles
    rigidBodyState.C

\*---------------------------------------------------------------------------*/

#ifndef rigidBodyState_H
#define rigidBodyState_H

#include "fvMeshFunctionObject.H"
#include "logFiles.H"
#include "rigidBodyMotion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class rigidBodyState Declaration
\*---------------------------------------------------------------------------*/

class rigidBodyState
:
    public fvMeshFunctionObject,
    public logFiles
{
    // Private Data

        word angleFormat_;

        //- List of the names of the rigid bodies
        wordList names_;


protected:

    // Protected Member Functions

        //- overloaded writeFileHeader from writeFile
        virtual void writeFileHeader(const label i = 0);


public:

    //- Runtime type information
    TypeName("rigidBodyState");


    // Constructors

        //- Construct from Time and dictionary
        rigidBodyState
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        rigidBodyState(const rigidBodyState&) = delete;


    //- Destructor
    virtual ~rigidBodyState();


    // Member Functions

        const RBD::rigidBodyMotion& motion() const;

        //- Read the rigidBodyState data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Write the rigidBodyState
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const rigidBodyState&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
