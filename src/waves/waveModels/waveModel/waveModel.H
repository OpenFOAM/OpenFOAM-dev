/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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
    Foam::waveModel

Description
    Generic base class for waves. Derived classes must implement field
    functions which return the elevation above the wave surface and the
    velocity field, both as a function of position.

SourceFiles
    waveModel.C

\*---------------------------------------------------------------------------*/

#ifndef waveModel_H
#define waveModel_H

#include "objectRegistry.H"
#include "dictionary.H"
#include "Function1.H"
#include "runTimeSelectionTables.H"
#include "vectorField.H"
#include "vector2DField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class waveModel Declaration
\*---------------------------------------------------------------------------*/

class waveModel
{
    // Private data

        //- Reference to the database
        const objectRegistry& db_;

        //- Peak-to-mean amplitude [m]
        autoPtr<Function1<scalar>> amplitude_;


public:

    //- Runtime type information
    TypeName("waveModel");


    // Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        waveModel,
        objectRegistry,
        (const objectRegistry& db, const dictionary& dict),
        (db, dict)
    );


    // Constructors

        //- Construct a copy
        waveModel(const waveModel& wave);

        //- Construct from a database and a dictionary
        waveModel(const objectRegistry& db, const dictionary& dict);

        //- Construct a clone
        virtual autoPtr<waveModel> clone() const = 0;


    // Selectors

        //- Select
        static autoPtr<waveModel> New
        (
            const objectRegistry& db,
            const dictionary& dict
        );

        //- Select
        static autoPtr<waveModel> New
        (
            const word& type,
            const objectRegistry& db,
            const dictionary& dict
        );


    //- Destructor
    virtual ~waveModel();


    // Member Functions

        // Access

            //- Get the amplitude
            scalar amplitude(const scalar t) const
            {
                return amplitude_->value(t);
            }

        //- Get the (scalar) value of gravity.
        scalar g() const;

        //- Get the wave elevation at a given time, mean velocity and local
        //  coordinates. Local x is aligned with the mean velocity.
        virtual tmp<scalarField> elevation
        (
            const scalar t,
            const scalar u,
            const scalarField& x
        ) const = 0;

        //- Get the wave velocity at a given time, mean velocity and local
        //  coordinates. Local x is aligned with the mean velocity, and z with
        //  negative gravity.
        virtual tmp<vector2DField> velocity
        (
            const scalar t,
            const scalar u,
            const vector2DField& xz
        ) const = 0;

        //- Get the wave pressure at a given time, mean velocity and local
        //  coordinates. Local x is aligned with the mean velocity, and z with
        //  negative gravity.
        virtual tmp<scalarField> pressure
        (
            const scalar t,
            const scalar u,
            const vector2DField& xz
        ) const = 0;

        //- Write
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
