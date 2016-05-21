/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2016 OpenFOAM Foundation
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
    Foam::functionObjects::vorticity

Group
    grpUtilitiesFunctionObjects

Description
    This function object calculates the vorticity, the curl of the velocity.

SeeAlso
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    vorticity.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_vorticity_H
#define functionObjects_vorticity_H

#include "fvMeshFunctionObject.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class vorticity Declaration
\*---------------------------------------------------------------------------*/

class vorticity
:
    public fvMeshFunctionObject
{
    // Private data

        //- Name of velocity field, default is "U"
        word UName_;

        //- Name of vorticity field
        word outputName_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        vorticity(const vorticity&);

        //- Disallow default bitwise assignment
        void operator=(const vorticity&);


public:

    //- Runtime type information
    TypeName("vorticity");


    // Constructors

        //- Construct from Time and dictionary
        vorticity
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~vorticity();


    // Member Functions

        //- Read the vorticity data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute(const bool postProcess = false);

        //- Calculate the vorticity and write
        virtual bool write(const bool postProcess = false);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
