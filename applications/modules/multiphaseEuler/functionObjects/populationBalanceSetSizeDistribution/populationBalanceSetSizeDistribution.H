/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::functionObjects::populationBalanceSetSizeDistribution

Description
    Sets the population balance size distribution by overwriting the values in
    the size-group fraction fields with values obtained by integrating a given
    distribution.

    This function will only execute in post-processing mode.

    Example of function object specification:
    \verbatim
    populationBalanceSetSizeDistribution
    {
        type        populationBalanceSetSizeDistribution;
        libs        ("libmultiphaseEulerFunctionObjects.so");

        populationBalance bubbles;

        //phase       air1; // <-- alternative to specifying
                            //     'populationBalance' if only a single phase's
                            //     size group fractions are to be set

        distribution
        {
            type        tabulatedDensity;
            Q           3;
            file        "constant/initialDistribution.bubbles";
        }
    }
    \endverbatim

Usage
    \table
        Property          | Description             | Required        | Default
        populationBalance | population balance name | if phase is not \
                                                      specified       |
        phase             | phase name              | if populationBalance is \
                                                      not specified   |
        distribution      | the distribution        | yes             |
    \endtable

See also
    Foam::diameterModels::populationBalanceModel
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObject

SourceFiles
    populationBalanceSetSizeDistribution.C

\*---------------------------------------------------------------------------*/

#ifndef populationBalanceSetSizeDistribution_H
#define populationBalanceSetSizeDistribution_H

#include "fvMeshFunctionObject.H"
#include "populationBalanceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
             Class populationBalanceSetSizeDistribution Declaration
\*---------------------------------------------------------------------------*/

class populationBalanceSetSizeDistribution
:
    public fvMeshFunctionObject
{
private:

    // Private Data

        //- Name of the population balance
        const word popBalName_;

        //- Name of the phase
        const word phaseName_;

        //- Distribution
        autoPtr<distribution> distribution_;


public:

    //- Runtime type information
    TypeName("populationBalanceSetSizeDistribution");


    // Constructors

        //- Construct from Time and dictionary
        populationBalanceSetSizeDistribution
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        populationBalanceSetSizeDistribution
        (
            const populationBalanceSetSizeDistribution&
        ) = delete;


    //- Destructor
    virtual ~populationBalanceSetSizeDistribution();


    // Member Functions

        //- Return the list of fields required
        virtual wordList fields() const
        {
            return wordList::null();
        }

        //- Execute. Does nothing.
        virtual bool execute();

        //- Calculate and write the size group fraction fields
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const populationBalanceSetSizeDistribution&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
