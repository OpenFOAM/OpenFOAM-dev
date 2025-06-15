/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024-2025 OpenFOAM Foundation
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
    Foam::functionObjects::reactionRates

Description
    Writes volume averaged reaction rates in [kmol/m^3/s] for
    each reaction into the file \<timeDir\>/reactionRates.dat

    Example of function object specification:
    \verbatim
    reactionRates1
    {
        type            reactionRates;
        libs            ("libchemistryModel.so");

        cellZone        all; // Or points, cellSet, cellZone

        //phase         <phaseName>; // Optional name of the phase

        //writeFields   false; // Optionally also write the rate fields. Note
                               // that this can create a lot of data.
    }
    \endverbatim

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::fvCellZone
    Foam::functionObjects::logFiles

SourceFiles
    reactionRates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_reactionRates_H
#define functionObjects_reactionRates_H

#include "fvMeshFunctionObject.H"
#include "fvCellZone.H"
#include "logFiles.H"
#include "basicChemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                     Class reactionRates Declaration
\*---------------------------------------------------------------------------*/

class reactionRates
:
    public fvMeshFunctionObject,
    public fvCellZone,
    public logFiles
{
    // Private Member Data

        //- Name of the phase
        const word phaseName_;

        //- Reference to the chemistry model
        const basicChemistryModel& chemistryModel_;

        //- Whether or not to write the reaction rate fields
        bool writeFields_;


    // Private Member Functions

        //- File header information
        virtual void writeFileHeader(const label i);


public:

    //- Runtime type information
    TypeName("reactionRates");


    // Constructors

        //- Construct from Time and dictionary
        reactionRates
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        reactionRates(const reactionRates&) = delete;


    //- Destructor
    virtual ~reactionRates();


    // Member Functions

        //- Read the reactionRates data
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const
        {
            return wordList::null();
        }

        //- Do nothing
        virtual bool execute();

        //- Write the reaction rates
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const reactionRates&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
