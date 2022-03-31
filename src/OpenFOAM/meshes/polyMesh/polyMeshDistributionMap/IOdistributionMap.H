/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2022 OpenFOAM Foundation
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
    Foam::IOdistributionMap

Description
    IOdistributionMap is derived from distributionMap and
    IOobject to give the distributionMap
    automatic IO functionality via the objectRegistry.

SourceFiles
    IOdistributionMap.C

\*---------------------------------------------------------------------------*/

#ifndef IOdistributionMap_H
#define IOdistributionMap_H

#include "distributionMap.H"
#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class IOdistributionMap Declaration
\*---------------------------------------------------------------------------*/

class IOdistributionMap
:
    public regIOobject,
    public distributionMap
{

public:

    //- Runtime type information
    TypeName("distributionMap");

    // Constructors

        //- Construct given an IOobject
        IOdistributionMap(const IOobject&);

        //- Construct given an IOobject and distributionMap
        IOdistributionMap(const IOobject&, const distributionMap&);

        //- Move constructor transferring the distributionMap contents
        IOdistributionMap(const IOobject&, distributionMap&&);


    //- Destructor
    virtual ~IOdistributionMap();


    // Member Functions

        //- ReadData function required for regIOobject read operation
        virtual bool readData(Istream&);

        //- WriteData function required for regIOobject write operation
        virtual bool writeData(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
