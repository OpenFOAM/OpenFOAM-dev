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
    Foam::consumptionSpeed

Description
    Correlation function for laminar consumption speed obtained from flamelet
    solution at increasing strain rates.

SourceFiles
    consumptionSpeed.C

\*---------------------------------------------------------------------------*/

#ifndef consumptionSpeed_H
#define consumptionSpeed_H

#include "IOdictionary.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class consumptionSpeed Declaration
\*---------------------------------------------------------------------------*/

class consumptionSpeed
{
    // Private Data

        //- Maximum consumption speed
        scalar omega0_;

        //- Exponential factor
        scalar eta_;

        //- Extinction strain
        scalar sigmaExt_;

        //- Minimum consumption speed
        scalar omegaMin_;


    // Private Member Functions

        //- Return consumption rate
        scalar omega0Sigma(scalar sigma, scalar a) const;


public:

    //- Runtime type information
    TypeName("consumptionSpeed");


    // Constructors

        //- Construct from dictionary
        consumptionSpeed(const dictionary& dict);

        //- Disallow default bitwise copy construction
        consumptionSpeed(const consumptionSpeed&) = delete;


    //- Destructor
    virtual ~consumptionSpeed();


    // Member Functions

        //- Return speed consumption rate temp
        tmp<volScalarField> omega0Sigma(const volScalarField& sigma);


        // Access functions

            scalar omega0() const
            {
                return omega0_;
            }

            scalar eta() const
            {
                return eta_;
            }

            scalar sigmaExt() const
            {
                return sigmaExt_;
            }

            scalar omegaMin() const
            {
                return omegaMin_;
            }


    // IO

        //- Update properties
        void read(const dictionary& dict);


    // Member Operators

         //- Disallow default bitwise assignment
        void operator=(const consumptionSpeed&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
