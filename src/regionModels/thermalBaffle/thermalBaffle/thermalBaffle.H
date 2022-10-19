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
    Foam::regionModels::thermalBaffle

Description
    Thermal baffle region model

SourceFiles
    thermalBaffle.C
    thermalBaffleI.H

\*---------------------------------------------------------------------------*/

#ifndef thermalBaffle_H
#define thermalBaffle_H

#include "regionModel.H"
#include "solidThermo.H"
#include "solidThermophysicalTransportModel.H"
#include "radiationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{

/*---------------------------------------------------------------------------*\
                        Class thermalBaffle Declaration
\*---------------------------------------------------------------------------*/

class thermalBaffle
:
    public regionModel
{
    // Private Member Functions

        //- Calculate the baffle mesh thickness
        tmp<volScalarField::Internal> calcThickness() const;


protected:

    // Protected data

        // Geometric parameters

            //- Baffle thickness
            volScalarField::Internal thickness_;


        // Solution parameters

            //- Number of non orthogonal correctors
            label nNonOrthCorr_;


        // Thermo properties

            //- Solid thermo
            autoPtr<solidThermo> thermo_;


        // Thermophysical transport

            //- Pointer to the solid thermophysical transport model
            autoPtr<solidThermophysicalTransportModel> thermophysicalTransport_;


        // Source term fields

            //- Surface energy source [J/m2/s]
            volScalarField Qs_;

            //- Volumetric energy source [J/m3/s]
            volScalarField Q_;


        // Sub models

            //- Pointer to radiation model
            autoPtr<radiationModel> radiation_;


    // Protected member functions

        // IO

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
            const solidThermo& thermo() const;


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

} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
