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
    Foam::DsmcBaseCloud

Description
    Virtual abstract base class for templated DsmcCloud

SourceFiles
    DsmcBaseCloud.C

\*---------------------------------------------------------------------------*/

#ifndef DsmcBaseCloud_H
#define DsmcBaseCloud_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class DsmcBaseCloud Declaration
\*---------------------------------------------------------------------------*/

class DsmcBaseCloud
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        DsmcBaseCloud(const DsmcBaseCloud&);

        //- Disallow default bitwise assignment
        void operator=(const DsmcBaseCloud&);


public:

    //- Runtime type information
    TypeName("DsmcBaseCloud");

    // Constructors

        //- Null constructor
        DsmcBaseCloud();

    //- Destructor
    virtual ~DsmcBaseCloud();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
