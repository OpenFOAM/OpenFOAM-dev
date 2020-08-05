/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2020 OpenFOAM Foundation
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
    Writes the domain averaged reaction rates for each specie for each reaction
    into the file \<timeDir\>/specieReactionRates.dat

See also
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::volRegion
    Foam::functionObjects::logFiles

SourceFiles
    specieReactionRates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_specieReactionRates_H
#define functionObjects_specieReactionRates_H

#include "fvMeshFunctionObject.H"
#include "volRegion.H"
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
    public volRegion,
    public logFiles
{
    // Private Member Data

        //- Reference to the chemistry model
        const basicChemistryModel& chemistryModel_;


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
