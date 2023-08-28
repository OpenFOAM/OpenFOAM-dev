/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2023 OpenFOAM Foundation
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
    Foam::genericFieldBase

Description
    Base class for generic field types. Facilitates down-casting so that the
    actual type can be queried.

SourceFiles
    genericFieldBase.C

\*---------------------------------------------------------------------------*/

#ifndef genericFieldBase_H
#define genericFieldBase_H

#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class genericFieldBase Declaration
\*---------------------------------------------------------------------------*/

class genericFieldBase
{
    // Private Data

        //- Actual field type name
        word actualTypeName_;


public:

    //- Runtime type information
    TypeName("generic");


    // Constructors

        //- Construct null
        genericFieldBase();

        //- Construct from actual type name
        genericFieldBase(const word& actualTypeName);


    //- Destructor
    virtual ~genericFieldBase();


    // Member Functions

        //- Return the actual type name
        inline const word& actualTypeName() const
        {
            return actualTypeName_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
