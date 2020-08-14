/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::IOPosition

Description
    Helper IO class to read and write particle positions

SourceFiles
    IOPosition.C

\*---------------------------------------------------------------------------*/

#ifndef IOPosition_H
#define IOPosition_H

#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class IOPositionName Declaration
\*---------------------------------------------------------------------------*/

TemplateName(IOPosition);


/*---------------------------------------------------------------------------*\
                         Class IOPosition Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class IOPosition
:
    public regIOobject,
    public IOPositionName
{

    // Private Data

        //- Reference to the cloud
        const CloudType& cloud_;


public:

    //- Type information

        using IOPositionName::typeName;

        virtual const word& type() const
        {
            return IOPositionName::typeName;
        }


    // Constructors

        //- Construct from cloud
        IOPosition(const CloudType&);


    // Member Functions

        //- Inherit readData from regIOobject
        using regIOobject::readData;

        virtual void readData(Istream&, CloudType&);

        virtual bool write(const bool write = true) const;

        virtual bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "IOPosition.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
