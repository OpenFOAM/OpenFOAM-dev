/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2022 OpenFOAM Foundation
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
    Model to describe the dependence of saturation pressure on temperature, and
    vice versa

SourceFiles
    saturationModel.C
    saturationModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef saturationModel_H
#define saturationModel_H

#include "volFields.H"
#include "dictionary.H"
#include "phaseInterface.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class saturationModel Declaration
\*---------------------------------------------------------------------------*/

class saturationModel
:
    public regIOobject
{
    // Private Data

        //- The interface
        const phaseInterface interface_;


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
            const dictionary& dict,
            const phaseInterface& interface
        ),
        (dict, interface)
    );


    // Constructors

        //- Construct from a dictionary and an interface
        saturationModel
        (
            const dictionary& dict,
            const phaseInterface& interface
        );


    // Selectors

        static autoPtr<saturationModel> New
        (
            const dictionary& dict,
            const phaseInterface& interface,
            const bool outer=true
        );


    //- Destructor
    virtual ~saturationModel();


    // Member Functions

        //- Access the interface
        const phaseInterface& interface() const;

        //- Saturation pressure
        virtual tmp<volScalarField> pSat
        (
            const volScalarField& T
        ) const = 0;

        //- Saturation pressure derivative w.r.t. temperature
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

        //- Dummy write for regIOobject
        bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
