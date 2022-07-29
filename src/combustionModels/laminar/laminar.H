/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2022 OpenFOAM Foundation
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
    Foam::combustionModels::laminar

Description
    Laminar combustion model.

SourceFiles
    laminar.C

\*---------------------------------------------------------------------------*/

#ifndef laminar_H
#define laminar_H

#include "autoPtr.H"
#include "combustionModel.H"
#include "basicChemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                           Class laminar Declaration
\*---------------------------------------------------------------------------*/

class laminar
:
    public combustionModel
{
    // Private Data

        //- Run chemistry correction on every outer iteration. Default false.
        bool outerCorrect_;

        //- Integrate reaction rate over the time-step
        //  using the selected ODE solver. Default true.
        bool integrateReactionRate_;

        //- Maximum integration time permitted. Default vGreat.
        scalar maxIntegrationTime_;

        //- The time index of the last correction
        label timeIndex_;

protected:

    // Protected Data

        //- Pointer to chemistry model
        autoPtr<basicChemistryModel> chemistryPtr_;


public:

    //- Runtime type information
    TypeName("laminar");


    // Constructors

        //- Construct from components
        laminar
        (
            const word& modelType,
            const fluidMulticomponentThermo& thermo,
            const compressibleMomentumTransportModel& turb,
            const word& combustionProperties
        );

        //- Disallow default bitwise copy construction
        laminar(const laminar&);


    //- Destructor
    virtual ~laminar();


    // Member Functions

        //- Correct combustion rate
        virtual void correct();

        //- Fuel consumption rate matrix.
        virtual tmp<fvScalarMatrix> R(volScalarField& Y) const;

        //- Heat release rate [kg/m/s^3]
        virtual tmp<volScalarField> Qdot() const;

        //- Update properties from given dictionary
        virtual bool read();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const laminar&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
