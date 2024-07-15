/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::SuModel

Description
    Base-class for all Su models used by the b-Xi combustion model

    References:
    \verbatim
        Weller, H. G. (1993).
        The development of a new flame area combustion model
        using conditional averaging.
        Thermo-fluids section report TF 9307.

        Weller, H. G., Tabor, G., Gosman, A. D., & Fureby, C. (1998, January).
        Application of a flame-wrinkling LES combustion model
        to a turbulent mixing layer.
        In Symposium (International) on combustion
        (Vol. 27, No. 1, pp. 899-907). Elsevier.
    \endverbatim

SourceFiles
    SuModel.C

\*---------------------------------------------------------------------------*/

#ifndef SuModel_H
#define SuModel_H

#include "psiuMulticomponentThermo.H"
#include "fluidThermoThermophysicalTransportModels.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class SuModel Declaration
\*---------------------------------------------------------------------------*/

class SuModel
{

protected:

    // Protected data

        //- Thermo
        const psiuMulticomponentThermo& thermo_;

        //- Thermo-physical transport
        const fluidThermoThermophysicalTransportModel& thermoTransport_;

        //- Turbulence
        const compressibleMomentumTransportModel& turbulence_;

        //- Laminar flame speed
        volScalarField Su_;

        //- Update coefficients from given dictionary
        virtual bool readCoeffs(const dictionary& dict) = 0;


public:

    //- Runtime type information
    TypeName("SuModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            SuModel,
            dictionary,
            (
                const dictionary& dict,
                const psiuMulticomponentThermo& thermo,
                const fluidThermoThermophysicalTransportModel& turbulence
            ),
            (
                dict,
                thermo,
                turbulence
            )
        );


    // Constructors

        //- Construct from components
        SuModel
        (
            const psiuMulticomponentThermo& thermo,
            const fluidThermoThermophysicalTransportModel& turbulence
        );

        //- Disallow default bitwise copy construction
        SuModel(const SuModel&) = delete;


    // Selectors

        //- Return a reference to the selected Su model
        static autoPtr<SuModel> New
        (
            const dictionary& combustionProperties,
            const psiuMulticomponentThermo& thermo,
            const fluidThermoThermophysicalTransportModel& turbulence
        );


    //- Destructor
    virtual ~SuModel();


    // Member Functions

        //- Return the laminar flame speed
        virtual const volScalarField& Su() const
        {
            return Su_;
        }

        //- Correct the laminar flame speed
        virtual void correct() = 0;

        //- Update properties from the given combustionProperties dictionary
        bool read(const dictionary& combustionProperties);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const SuModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
