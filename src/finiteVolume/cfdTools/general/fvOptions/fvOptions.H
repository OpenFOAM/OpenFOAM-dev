/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::fv::options

Description
    Finite-volume options

SourceFiles
    options.C

\*---------------------------------------------------------------------------*/

#ifndef fvOptions_H
#define fvOptions_H

#include "fvOptionList.H"
#include "IOdictionary.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class options Declaration
\*---------------------------------------------------------------------------*/

class options
:
    public IOdictionary,
    public optionList
{
    // Private Member Functions

        //- Create IO object if dictionary is present
        IOobject createIOobject(const fvMesh& mesh) const;

        //- Disallow default bitwise copy construct
        options(const options&);

        //- Disallow default bitwise assignment
        void operator=(const options&);


public:

    // Declare name of the class and its debug switch
    ClassName("fvOptions");


    // Constructors

        //- Construct from components with list of field names
        options(const fvMesh& mesh);

        //- Construct fvOptions and register to datbase if not present
        //  otherwise lookup and return
        static options& New(const fvMesh& mesh);


    //- Destructor
    virtual ~options()
    {}


    // Member Functions

        //- Inherit read from optionList
        using optionList::read;

        //- Read dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
