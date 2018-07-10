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
    Foam::pairPotentials::maitlandSmith

Description

    Reference:
    \verbatim
        Maitland, G. C., & Smith, E. B. (1973).
        A simplified representation of intermolecular potential energy.
        Chemical Physics Letters, 22(3), 443-446.
    \endverbatim

    Parameters for other monoatomics from:
    \verbatim
        Maitland, G. C., Rigby, M., Smith, E., Wakeham, W. (1981).
        Intermolecular forces: Their origin and determination.
        Oxford: Clarendon Press.
    \endverbatim

SourceFiles
    maitlandSmith.C

\*---------------------------------------------------------------------------*/

#ifndef maitlandSmith_H
#define maitlandSmith_H

#include "pairPotential.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace pairPotentials
{

/*---------------------------------------------------------------------------*\
                       Class maitlandSmith Declaration
\*---------------------------------------------------------------------------*/

class maitlandSmith
:
    public pairPotential
{
    // Private data

        dictionary maitlandSmithCoeffs_;

        scalar m_;
        scalar gamma_;
        scalar rm_;
        scalar epsilon_;


public:

    //- Runtime type information
    TypeName("maitlandSmith");


    // Constructors

        //- Construct from components
        maitlandSmith
        (
            const word& name,
            const dictionary& pairPotentialProperties
        );


    //- Destructor
    ~maitlandSmith()
    {}


    // Member Functions

        scalar unscaledEnergy(const scalar r) const;

        //- Read dictionary
        bool read(const dictionary& pairPotentialProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pairPotentials
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
