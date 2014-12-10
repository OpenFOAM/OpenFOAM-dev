/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::reactingMultiphaseCloud

Description
    Virtual abstract base class for templated reactingMultiphaseCloud

SourceFiles
    reactingMultiphaseCloud.C

\*---------------------------------------------------------------------------*/

#ifndef reactingMultiphaseCloud_H
#define reactingMultiphaseCloud_H

#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class reactingMultiphaseCloud Declaration
\*---------------------------------------------------------------------------*/

class reactingMultiphaseCloud
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        reactingMultiphaseCloud(const reactingMultiphaseCloud&);

        //- Disallow default bitwise assignment
        void operator=(const reactingMultiphaseCloud&);


public:

    //- Runtime type information
    TypeName("reactingMultiphaseCloud");

    // Constructors

        //- Null constructor
        reactingMultiphaseCloud();


    //- Destructor
    virtual ~reactingMultiphaseCloud();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
