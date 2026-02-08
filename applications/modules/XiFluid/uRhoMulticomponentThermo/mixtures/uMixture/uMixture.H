/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::uMixture

Description
    Unburnt gas mixture base class

SourceFiles
    uMixture.C

\*---------------------------------------------------------------------------*/

#ifndef uMixture_H
#define uMixture_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class uMixture Declaration
\*---------------------------------------------------------------------------*/

class uMixture
{
public:

    // Constructors

        //- Construct
        uMixture()
        {}

        //- Disallow default bitwise copy construction
        uMixture(const uMixture&) = delete;


    //- Destructor
    virtual ~uMixture()
    {}


    // Member Functions

        //- Return the mixture type name
        virtual const word& mixtureType() const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const uMixture&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
