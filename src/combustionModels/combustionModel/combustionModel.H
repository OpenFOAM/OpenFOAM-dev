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
    Foam::combustionModel

Description
    Base class for combustion models

SourceFiles
    combustionModel.C

\*---------------------------------------------------------------------------*/

#ifndef combustionModel_H
#define combustionModel_H

#include "IOdictionary.H"
#include "compressibleMomentumTransportModels.H"
#include "fluidMulticomponentThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class combustionModel Declaration
\*---------------------------------------------------------------------------*/

class combustionModel
:
    public IOdictionary
{
    // Private Member Functions

        //- Construct the base IO object
        IOobject createIOobject
        (
            const fluidMulticomponentThermo& thermo,
            const word& combustionProperties
        ) const;


protected:

    // Protected data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Reference to the thermo
        const fluidMulticomponentThermo& thermo_;

        //- Reference to the turbulence model
        const compressibleMomentumTransportModel& turb_;

        //- Dictionary of the model
        dictionary coeffs_;

        //- Model type
        const word modelType_;


public:

    //- Runtime type information
    TypeName("combustionModel");

    //- Default combustionProperties dictionary name
    static const word combustionPropertiesName;


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        combustionModel,
        dictionary,
        (
            const word& modelType,
            const fluidMulticomponentThermo& thermo,
            const compressibleMomentumTransportModel& turb,
            const word& combustionProperties
        ),
        (modelType, thermo, turb, combustionProperties)
    );


    // Constructors

        //- Construct from components
        combustionModel
        (
            const word& modelType,
            const fluidMulticomponentThermo& thermo,
            const compressibleMomentumTransportModel& turb,
            const word& combustionProperties = combustionPropertiesName
        );

        //- Disallow default bitwise copy construction
        combustionModel(const combustionModel&) = delete;


    // Selectors

        //- Select using thermo and turbulence
        static autoPtr<combustionModel> New
        (
            const fluidMulticomponentThermo& thermo,
            const compressibleMomentumTransportModel& turb,
            const word& combustionProperties = combustionPropertiesName
        );


    //- Destructor
    virtual ~combustionModel();


    // Member Functions

        //- Return const access to the mesh database
        inline const fvMesh& mesh() const;

        //- Return const access to the thermo
        inline const fluidMulticomponentThermo& thermo() const;

        //- Return access to turbulence
        inline const compressibleMomentumTransportModel& turbulence() const;

        //- Return const access to rho
        inline const volScalarField& rho() const;

        //- Return const access to phi
        inline tmp<surfaceScalarField> phi() const;

        //- Return const dictionary of the model
        inline const dictionary& coeffs() const;

        //- Correct combustion rate
        virtual void correct() = 0;

        //- Fuel consumption rate matrix, i.e. source term for fuel equation
        virtual tmp<fvScalarMatrix> R(volScalarField& Y) const = 0;

        //- Heat release rate [kg/m/s^3]
        virtual tmp<volScalarField> Qdot() const = 0;

        //- Update properties from given dictionary
        virtual bool read();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const combustionModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "combustionModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
