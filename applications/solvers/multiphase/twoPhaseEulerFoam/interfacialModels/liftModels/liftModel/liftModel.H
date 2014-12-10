/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 OpenFOAM Foundation
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
    Foam::liftModel

Description

SourceFiles
    liftModel.C
    newLiftModel.C

\*---------------------------------------------------------------------------*/

#ifndef liftModel_H
#define liftModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class phasePair;

/*---------------------------------------------------------------------------*\
                           Class liftModel Declaration
\*---------------------------------------------------------------------------*/

class liftModel
{
protected:

    // Protected data

        //- Phase pair
        const phasePair& pair_;


public:

    //- Runtime type information
    TypeName("liftModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            liftModel,
            dictionary,
            (
                const dictionary& dict,
                const phasePair& pair
            ),
            (dict, pair)
        );


    // Static data members

        //- Force dimensions
        static const dimensionSet dimF;


    // Constructors

        //- Construct from a dictionary and a phase pair
        liftModel
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~liftModel();


    // Selectors

        static autoPtr<liftModel> New
        (
            const dictionary& dict,
            const phasePair& pair
        );


    // Member Functions

        //- Lift coefficient
        virtual tmp<volScalarField> Cl() const = 0;

        //- Lift force
        virtual tmp<volVectorField> F() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
