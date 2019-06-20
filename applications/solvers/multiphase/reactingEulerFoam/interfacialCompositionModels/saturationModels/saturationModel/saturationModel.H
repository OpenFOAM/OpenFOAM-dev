/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2019 OpenFOAM Foundation
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
    Foam::saturationModel

Description

SourceFiles
    saturationModel.C
    newSaturationModel.C

\*---------------------------------------------------------------------------*/

#ifndef saturationModel_H
#define saturationModel_H

#include "IOdictionary.H"
#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class saturationModel Declaration
\*---------------------------------------------------------------------------*/

class saturationModel
:
    public IOdictionary
{
public:

    //- Runtime type information
    TypeName("saturationModel");


    //- Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        saturationModel,
        dictionary,
        (
            const dictionary& dict, const objectRegistry& db
        ),
        (dict, db)
    );


    // Constructors

        //- Construct null
        saturationModel(const objectRegistry& db);

        //- Disallow default bitwise copy construction
        saturationModel(const saturationModel&) = delete;


    // Selectors

        //- Select null constructed
        static autoPtr<saturationModel> New
        (
            const dictionary& dict,
            const objectRegistry& db
        );


    //- Destructor
    virtual ~saturationModel();


    // Member Functions

        //- Saturation pressure
        virtual tmp<volScalarField> pSat
        (
            const volScalarField& T
        ) const = 0;

        //- Saturation pressure derivetive w.r.t. temperature
        virtual tmp<volScalarField> pSatPrime
        (
            const volScalarField& T
        ) const = 0;

        //- Natural log of the saturation pressure
        virtual tmp<volScalarField> lnPSat
        (
            const volScalarField& T
        ) const = 0;

        //- Saturation temperature
        virtual tmp<volScalarField> Tsat
        (
            const volScalarField& p
        ) const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const saturationModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
