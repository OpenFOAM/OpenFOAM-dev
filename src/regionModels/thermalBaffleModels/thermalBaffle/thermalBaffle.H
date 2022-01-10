/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::regionModels::thermalBaffleModels::thermalBaffle

Description
    2D thermal baffle

SourceFiles
    thermalBaffle.C
    thermalBaffleI.H

\*---------------------------------------------------------------------------*/

#ifndef thermalBaffle_H
#define thermalBaffle_H

#include "thermalBaffleModel.H"
#include "volFieldsFwd.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace thermalBaffleModels
{


/*---------------------------------------------------------------------------*\
                        Class thermalBaffle Declaration
\*---------------------------------------------------------------------------*/

class thermalBaffle
:
    public thermalBaffleModel
{
    // Private Member Functions

        //- Initialise thermalBaffle
        void init();


protected:

    // Protected data

        // Solution parameters

            //- Number of non orthogonal correctors
            label nNonOrthCorr_;


        // Thermo properties

            //- Solid thermo
            autoPtr<solidThermo> thermo_;

            //- Enthalpy/internal energy
            volScalarField& he_;


        // Source term fields

            //- Surface energy source  / [J/m2/s]
            volScalarField Qs_;

            //- Volumetric energy source  / [J/m3/s]
            volScalarField Q_;


        // Sub models

            //- Pointer to radiation model
            autoPtr<radiationModel> radiation_;


    // Protected member functions

        //- Read control parameters IO dictionary
        virtual bool read();

        //- Read control parameters from dictionary
        virtual bool read(const dictionary& dict);


        // Equations

            //- Solve energy equation
            void solveEnergy();


public:

    //- Runtime type information
    TypeName("thermalBaffle");


    // Constructors

        //- Construct from components
        thermalBaffle(const word& modelType, const fvMesh& mesh);

        //- Construct from components and dict
        thermalBaffle
        (
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        thermalBaffle(const thermalBaffle&) = delete;


    //- Destructor
    virtual ~thermalBaffle();


    // Member Functions

        // Thermo properties

            //- Return const reference to the solidThermo
            virtual const solidThermo& thermo() const;


            // Fields

                //- Return the film specific heat capacity [J/kg/K]
                virtual const tmp<volScalarField> Cp() const;

                //- Return solid absorptivity [1/m]
                virtual const volScalarField& kappaRad() const;

                //- Return temperature [K]
                virtual const volScalarField& T() const;

                //- Return density [Kg/m^3]
                virtual const volScalarField& rho() const;

                //- Return thermal conductivity [W/m/K]
                virtual const volScalarField& kappa() const;


            // Helper functions

                //- Return sensible enthalpy/internal energy
                //  as a function of temperature
                //  for a patch
                inline tmp<scalarField> he
                (
                    const scalarField& p,
                    const scalarField& T,
                    const label patchi
                ) const;

                //- Return sensible enthalpy/internal energy
                inline tmp<volScalarField> he() const;


        // Evolution

            //- Pre-evolve  thermal baffle
            virtual void preEvolveRegion();

            //- Evolve the thermal baffle
            virtual void evolveRegion();


       // I-O

            //- Provide some feedback
            virtual void info();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const thermalBaffle&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace thermalBaffleModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "thermalBaffleI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
