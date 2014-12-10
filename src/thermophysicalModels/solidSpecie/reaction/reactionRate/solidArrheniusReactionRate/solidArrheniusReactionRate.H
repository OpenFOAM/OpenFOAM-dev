/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::solidArrheniusReactionRate

Description
    Arrhenius reaction rate for solids

SourceFiles
    solidArrheniusReactionRateI.H

\*---------------------------------------------------------------------------*/

#ifndef solidArrheniusReactionRate_H
#define solidArrheniusReactionRate_H

#include "scalarField.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class solidArrheniusReactionRate Declaration
\*---------------------------------------------------------------------------*/

class solidArrheniusReactionRate
{
    // Private data

        scalar A_;
        scalar Ta_;
        scalar Tcrit_;


public:

    // Constructors

        //- Construct from components
        inline solidArrheniusReactionRate
        (
            const scalar A,
            const scalar Ta,
            const scalar Tcrit
        );

        //- Construct from Istream
        inline solidArrheniusReactionRate
        (
            const speciesTable& species,
            Istream& is
        );

        //- Construct from dictionary
        inline solidArrheniusReactionRate
        (
            const speciesTable& species,
            const dictionary& dict
        );


    //- Destructor
    virtual ~solidArrheniusReactionRate()
    {}


    // Member Functions

        //- Return the type name
        static word type()
        {
            return "Arrhenius";
        }

        inline scalar operator()
        (
            const scalar p,
            const scalar T,
            const scalarField& c
        ) const;


         //- Write to stream
        inline void write(Ostream& os) const;


    // Ostream Operator

        inline friend Ostream& operator<<
        (
            Ostream&,
            const solidArrheniusReactionRate&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "solidArrheniusReactionRateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
