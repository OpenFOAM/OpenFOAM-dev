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
    Foam::wallBoilingModels::departureDiameterModel

Description
    Base class for bubble departure diameter models

SourceFiles
    departureDiameterModel.C
    departureDiameterModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef departureDiameterModel_H
#define departureDiameterModel_H

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

#include "phaseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{

/*---------------------------------------------------------------------------*\
                   Class departureDiameterModel Declaration
\*---------------------------------------------------------------------------*/

class departureDiameterModel
{
public:

    //- Runtime type information
    TypeName("departureDiameterModel");


    //- Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        departureDiameterModel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );


    // Constructors

        //- Construct null
        departureDiameterModel();

        //- Copy constructor
        departureDiameterModel(const departureDiameterModel&);

        //- Construct and return a clone
        virtual autoPtr<departureDiameterModel> clone() const = 0;


    // Selectors

        //- Select null constructed
        static autoPtr<departureDiameterModel> New(const dictionary& dict);


    //- Destructor
    virtual ~departureDiameterModel();


    // Member Functions

        //- Calculate and return the departure diameter field
        virtual tmp<scalarField> dDeparture
        (
            const phaseModel& liquid,
            const phaseModel& vapor,
            const label patchi,
            const scalarField& Tl,
            const scalarField& Tsatw,
            const scalarField& L
        ) const = 0;

        virtual void write(Ostream& os) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const departureDiameterModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace wallBoilingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
