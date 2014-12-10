/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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
    Foam::helpFunctionObject

Description
    This class provides help for functionObjects.  When no additional arguments
    are given, the utility outputs all known function objects.

    Usage:
    \plaintable
        \c -browse \<patch type\> | open documentation for object in browser
    \endplaintable

SeeAlso
    Foam::helpType
    Foam::functionObject

SourceFiles
    helpFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef helpFunctionObject_H
#define helpFunctionObject_H

#include "helpType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace helpTypes
{

/*---------------------------------------------------------------------------*\
                     Class helpFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class helpFunctionObject
:
    public helpType
{

public:

    //- Runtime type information
    TypeName("helpFunctionObject");

    //- Constructor
    helpFunctionObject();

    //- Destructor
    virtual ~helpFunctionObject();


    // Member Functions

        //- Initialise - typically setting static variables,
        //  e.g. command line arguments
        virtual void init();

        //- Execute the help
        virtual void execute(const argList& args, const fvMesh& mesh);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace helpTypes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
