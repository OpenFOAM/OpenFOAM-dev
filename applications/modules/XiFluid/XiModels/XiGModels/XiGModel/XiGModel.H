/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::XiGModel

Description
    Base class for flame wrinkling \c Xi generation rate models

    See Foam::XiModel for more details on flame wrinkling modelling.

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
    XiGModel.C

\*---------------------------------------------------------------------------*/

#ifndef XiGModel_H
#define XiGModel_H

#include "IOdictionary.H"
#include "psiuMulticomponentThermo.H"
#include "fluidThermoThermophysicalTransportModels.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class XiGModel Declaration
\*---------------------------------------------------------------------------*/

class XiGModel
{

protected:

    // Protected data

        //- Thermo
        const psiuMulticomponentThermo& thermo_;

        //- Thermo-physical transport
        const fluidThermoThermophysicalTransportModel& thermoTransport_;

        //- Turbulence
        const compressibleMomentumTransportModel& turbulence_;

        const volScalarField& Su_;

        //- Update coefficients from given dictionary
        virtual bool readCoeffs(const dictionary& dict) = 0;


public:

    //- Runtime type information
    TypeName("XiGModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            XiGModel,
            dictionary,
            (
                const dictionary& dict,
                const psiuMulticomponentThermo& thermo,
                const fluidThermoThermophysicalTransportModel& turbulence,
                const volScalarField& Su
            ),
            (
                dict,
                thermo,
                turbulence,
                Su
            )
        );


    // Constructors

        //- Construct from components
        XiGModel
        (
            const psiuMulticomponentThermo& thermo,
            const fluidThermoThermophysicalTransportModel& turbulence,
            const volScalarField& Su
        );

        //- Disallow default bitwise copy construction
        XiGModel(const XiGModel&) = delete;


    // Selectors

        //- Return a reference to the selected XiG model
        static autoPtr<XiGModel> New
        (
            const dictionary& XiProperties,
            const psiuMulticomponentThermo& thermo,
            const fluidThermoThermophysicalTransportModel& turbulence,
            const volScalarField& Su
        );


    //- Destructor
    virtual ~XiGModel();


    // Member Functions

        //- Return the flame-wrinkling generation rate
        virtual tmp<volScalarField> G() const = 0;

        //- Update properties from the given XiProperties dictionary
        bool read(const dictionary& XiProperties);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const XiGModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
