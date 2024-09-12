/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2024 OpenFOAM Foundation
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
    Foam::functionObjects::specieReactionRates

Description
    Writes volume averaged reaction rates in [kg/m^3/s] for each specie and
    each reaction into the file \<timeDir\>/specieReactionRates.dat

    Example of function object specification:
    \verbatim
    specieReactionRates1
    {
        type            specieReactionRates;
        libs            ("libchemistryModel.so");

        select          all; // Or points, cellSet, cellZone

        //phase         <phaseName>; // Optional name of the phase

        //writeFields   false; // Optionally also write the rate fields. Note
                               // that this can create a lot of data.
    }
    \endverbatim

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::fvCellSet
    Foam::functionObjects::logFiles

SourceFiles
    specieReactionRates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_specieReactionRates_H
#define functionObjects_specieReactionRates_H

#include "fvMeshFunctionObject.H"
#include "fvCellSet.H"
#include "logFiles.H"
#include "basicChemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                     Class specieReactionRates Declaration
\*---------------------------------------------------------------------------*/

class specieReactionRates
:
    public fvMeshFunctionObject,
    public fvCellSet,
    public logFiles
{
    // Private Member Data

        //- Name of the phase
        const word phaseName_;

        //- Reference to the chemistry model
        const basicChemistryModel& chemistryModel_;

        //- Whether or not to write the specie reaction rate fields
        bool writeFields_;


    // Private Member Functions

        //- File header information
        virtual void writeFileHeader(const label i);


public:

    //- Runtime type information
    TypeName("specieReactionRates");


    // Constructors

        //- Construct from Time and dictionary
        specieReactionRates
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        specieReactionRates(const specieReactionRates&) = delete;


    //- Destructor
    virtual ~specieReactionRates();


    // Member Functions

        //- Read the specieReactionRates data
        virtual bool read(const dictionary&);

        //- Return the list of fields required
        virtual wordList fields() const
        {
            return wordList::null();
        }

        //- Do nothing
        virtual bool execute();

        //- Write the specie reaction rates
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const specieReactionRates&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
