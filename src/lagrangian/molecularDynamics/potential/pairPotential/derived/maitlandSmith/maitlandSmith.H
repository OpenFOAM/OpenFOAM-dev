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
    Foam::pairPotentials::maitlandSmith

Description
    Foam::pairPotentials::maitlandSmith

    From:
    \verbatim
        @ARTICLE{MA_Maitland_Smith,
        author = {{Maitland}, G.~C. and {Smith}, E.~B.},
        title = {A simplified representation of
        intermolecular potential energy},
        journal = {Chemical Physics Letters},
        year = 1973,
        month = oct,
        volume = 22,
        pages = {443-446},
        adsurl = {http://adsabs.harvard.edu/abs/1973CPL....22..443M},
        adsnote = {Provided by the SAO/NASA Astrophysics Data System}
        }
    \endverbatim

    Parameters for other monoatomics from:
    \verbatim
        @BOOK{MD_Maitland_Rigby_Smith_Wakeham,
        AUTHOR =       {Geoffrey C. Maitland and Maurice Rigby and
        E. Brian Smith and William A. Wakeham},
        TITLE =        {Intermolecular Forces: Their Origin and Determination},
        PUBLISHER =    {Oxford University Press},
        YEAR =         {1981}
        }
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
