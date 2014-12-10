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
    Foam::pairPotentials::azizChen

Description
    Foam::pairPotentials::azizChen

    From:
    \verbatim
        \article{MA_Aziz_Chen,
        author = {R. A. Aziz and H. H. Chen},
        collaboration = {},
        title = {An accurate intermolecular potential for argon},
        publisher = {AIP},
        year = {1977},
        journal = {The Journal of Chemical Physics},
        volume = {67},
        number = {12},
        pages = {5719-5726},
        url = {http://link.aip.org/link/?JCP/67/5719/1},
        doi = {10.1063/1.434827}
        }
    \endverbatim

SourceFiles
    azizChen.C

\*---------------------------------------------------------------------------*/

#ifndef azizChen_H
#define azizChen_H

#include "pairPotential.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace pairPotentials
{

/*---------------------------------------------------------------------------*\
                           Class azizChen Declaration
\*---------------------------------------------------------------------------*/

class azizChen
:
    public pairPotential
{
    // Private data

        dictionary azizChenCoeffs_;

        scalar epsilon_;
        scalar rm_;
        scalar A_;
        scalar alpha_;
        scalar C6_;
        scalar C8_;
        scalar C10_;
        scalar D_;
        scalar gamma_;


public:

    //- Runtime type information
    TypeName("azizChen");


    // Constructors

        //- Construct from components
        azizChen
        (
            const word& name,
            const dictionary& pairPotentialProperties
        );


    //- Destructor
    ~azizChen()
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
