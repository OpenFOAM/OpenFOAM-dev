/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2019 OpenFOAM Foundation
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
    Foam::IOporosityModelList

Description
    List of porosity models with IO functionality

SourceFiles
    IOporosityModelList.C

\*---------------------------------------------------------------------------*/

#ifndef IOporosityModelList_H
#define IOporosityModelList_H

#include "IOdictionary.H"
#include "porosityModelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class IOporosityModelList Declaration
\*---------------------------------------------------------------------------*/

class IOporosityModelList
:
    public IOdictionary,
    public porosityModelList
{
    // Private Member Functions

        //- Create IO object if dictionary is present
        IOobject createIOobject(const fvMesh& mesh) const;


public:

    // Constructors

        //- Construct from mesh
        IOporosityModelList(const fvMesh& mesh);

        //- Disallow default bitwise copy construction
        IOporosityModelList(const IOporosityModelList&) = delete;

        //- Destructor
        virtual ~IOporosityModelList()
        {}


    // Member Functions

        //- Read dictionary
        virtual bool read();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const IOporosityModelList&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
