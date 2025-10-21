/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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
    Foam::XiEqModel

Description
    Base class for equilibrium flame wrinkling \c XiEq models

SourceFiles
    XiEqModel.C

\*---------------------------------------------------------------------------*/

#ifndef XiEqModel_H
#define XiEqModel_H

#include "IOdictionary.H"
#include "ubRhoThermo.H"
#include "compressibleMomentumTransportModels.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class XiEqModel Declaration
\*---------------------------------------------------------------------------*/

class XiEqModel
{

protected:

    // Protected data

        //- Thermo
        const ubRhoThermo& thermo_;

        //- Turbulence
        const compressibleMomentumTransportModel& momentumTransport_;

        //- Laminar burning velocity
        const volScalarField& Su_;

        //- Update coefficients from given dictionary
        virtual bool readCoeffs(const dictionary& dict) = 0;


public:

    //- Runtime type information
    TypeName("XiEqModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            XiEqModel,
            dictionary,
            (
                const dictionary& dict,
                const ubRhoThermo& thermo,
                const compressibleMomentumTransportModel& momentumTransport,
                const volScalarField& Su
            ),
            (
                dict,
                thermo,
                momentumTransport,
                Su
            )
        );


    // Constructors

        //- Construct from components
        XiEqModel
        (
            const ubRhoThermo& thermo,
            const compressibleMomentumTransportModel& momentumTransport,
            const volScalarField& Su
        );

        //- Disallow default bitwise copy construction
        XiEqModel(const XiEqModel&) = delete;


    // Selectors

        //- Return a reference to the selected XiEq model
        static autoPtr<XiEqModel> New
        (
            const dictionary& XiProperties,
            const ubRhoThermo& thermo,
            const compressibleMomentumTransportModel& turbulence,
            const volScalarField& Su
        );


    //- Destructor
    virtual ~XiEqModel();


    // Member Functions

        //- Return the flame-wrinkling XiEq
        virtual tmp<volScalarField> XiEq() const = 0;

        //- Return the flame diffusivity
        virtual tmp<volScalarField> Db() const
        {
            return volScalarField::New
            (
                "Db",
                momentumTransport_.rho()*momentumTransport_.nut()
            );
        }

        //- Update properties from the given XiProperties dictionary
        bool read(const dictionary& XiProperties);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const XiEqModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
