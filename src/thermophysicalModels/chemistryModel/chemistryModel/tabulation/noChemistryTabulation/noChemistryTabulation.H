/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2022 OpenFOAM Foundation
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
    Foam::chemistryTabulationMethods::none

Description
    A chemistryTabulationMethod which does nothing to allow tabulation to be
    switched-off.

SourceFiles
    noChemistryTabulation.C

\*---------------------------------------------------------------------------*/

#ifndef noChemistryTabulation_H
#define noChemistryTabulation_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "chemistryTabulationMethod.H"

namespace Foam
{
namespace chemistryTabulationMethods
{

/*---------------------------------------------------------------------------*\
                            Class none Declaration
\*---------------------------------------------------------------------------*/

class none
:
    public chemistryTabulationMethod
{
public:

    //- Runtime type information
    TypeName("none");

    // Constructors

        //- Construct from dictionary
        none
        (
            const dictionary& chemistryProperties,
            const odeChemistryModel& chemistry
        );

        //- Disallow default bitwise copy construction
        none(const none&) = delete;


    // Destructor
    virtual ~none();


    // Member Functions

        //- Return false as this tabulation method does not tabulate
        virtual bool tabulates()
        {
            return false;
        }

        //- Find the closest stored leaf of phiQ and store the result in
        //  RphiQ or return false.
        virtual bool retrieve
        (
            const Foam::scalarField& phiq,
            scalarField& Rphiq
        )
        {
            return false;
        }

        // Add information to the tabulation.This function can grow an
        // existing point or add a new leaf to the binary tree Input : phiq
        // the new composition to store Rphiq the mapping of the new
        // composition point
        virtual label add
        (
            const scalarField& phiq,
            const scalarField& Rphiq,
            const label nActive,
            const label li,
            const scalar deltaT
        )
        {
            return false;
        }

        virtual bool update()
        {
            return false;
        }

        virtual void reset()
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace chemistryTabulationMethods
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
