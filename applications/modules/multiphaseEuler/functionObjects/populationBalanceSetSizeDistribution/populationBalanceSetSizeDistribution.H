/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024-2026 OpenFOAM Foundation
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
    the group fraction fields with values obtained by integrating a given
    distribution.

    If run in post-processing mode, then it will set the fields and write them
    out. If executed at run-time then it will set the fields at the start of
    the run but not then write them out.

    Example of function object specification:
    \verbatim
    populationBalanceSetSizeDistribution
    {
        type        populationBalanceSetSizeDistribution;
        libs        ("libmultiphaseEulerFunctionObjects.so");

        populationBalance bubbles;

        //phase       air1; // <-- alternative to specifying
                            //     'populationBalance' if only a single phase's
                            //     group fractions are to be set

        distribution
        {
            type        tabulatedDensity;
            Q           3;
            file        "constant/initialDistribution.bubbles";
        }
    }
    \endverbatim

    Example of function object specification, using different distributions in
    different zones:
    \verbatim
    populationBalanceSetSizeDistribution
    {
        type        populationBalanceSetSizeDistribution;
        libs        ("libmultiphaseEulerFunctionObjects.so");

        populationBalance bubbles;

        default
        {
            type        multiFixedValue;
            Q           3;
            file        "constant/defaultDistribution";
        }

        zones
        {
            myCellZone0
            {
                type        RosinRammler;
                Q           0;
                min         0.2 [mm];
                max         2 [mm];
                d           1 [mm];
                n           2;
            }

            myCellZone1
            {
                type        tabulatedDensity;
                Q           3;
                distribution
                (
                    (10e-6 0.0025)
                    (15e-6 0.0528)
                    (20e-6 0.2795)
                    ...
                );
            }
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
    Foam::populationBalanceModel
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

        //- Default distribution
        autoPtr<distribution> defaultDistribution_;

        //- The indices of the zones for which non-default distributions apply
        labelList zoneIndices_;

        //- The distributions applied to the zones
        PtrList<distribution> zoneDistributions_;


    // Private Member Functions

        //- Set a population balance field
        template<class EtaVArg>
        void setField
        (
            const populationBalanceModel& popBal,
            const EtaVArg& etaVArg,
            volScalarField& f
        );

        //- Set the population balance fields
        void setFields();

        //- Write the population balance fields
        void writeFields();


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

        //- Calculate and write the group fraction fields
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
