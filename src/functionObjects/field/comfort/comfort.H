/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2020 OpenFOAM Foundation
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
    Foam::functionObjects::comfort

Description
    Calculates the thermal comfort quantities predicted mean vote (PMV) and
    predicted percentage of dissatisfaction (PPD) based on DIN ISO EN 7730:2005.

Usage
    \table
        Property      | Description                  | Required  | Default value
        clothing      | The insulation value of the cloth | no   | 0
        metabolicRate | The metabolic rate      | no        | 0.8
        extWork       | The external work        | no        | 0
        Trad          | Radiation temperature | no | -1
        relHumidity   | Relative humidity of the air | no | 50
        pSat          | Saturation pressure of water | no | -1
        tolerance     | Residual control for the cloth temperature | no | 1e-5
        maxClothIter  | Maximum number of iterations | no       | 0
        meanVelocity  | Use a constant mean velocity in the whole domain | no |\
        false
    \endtable

    \table
        Predicted Mean Vote (PMV)   | evaluation
        + 3                         | hot
        + 2                         | warm
        + 1                         | slightly warm
        + 0                         | neutral
        - 1                         | slightly cool
        - 2                         | cool
        - 3                         | cold
    \endtable

    \verbatim
    comfortAnalysis
    {
        type            comfort;
        libs            ("libfieldFunctionObjects.so");

        executeControl  writeTime;
        writeControl    writeTime;
    }
    \endverbatim

SourceFiles
    comfort.C

\*---------------------------------------------------------------------------*/

#ifndef comfort_H
#define comfort_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                           Class comfort Declaration
\*---------------------------------------------------------------------------*/

class comfort
:
    public fvMeshFunctionObject
{
    // Private Data

        //- Clothing [-]
        dimensionedScalar clothing_;

        //- Metabolic rate [kg/s^3]
        dimensionedScalar metabolicRate_;

        //- External work [kg/s^3]
        dimensionedScalar extWork_;

        //- Switch set to true if the radiation temperature is provided
        Switch TradSet_;

        //- Mean radiation temperature [K]
        dimensionedScalar Trad_;

        //- Relative humidity [percentage]
        dimensionedScalar relHumidity_;

        //- Saturation pressure of water [Pa]
        dimensionedScalar pSat_;

        //- Thermal insulation of clothing [W/m^2/K]
        dimensionedScalar Icl_;

        //- Prefactor of cloth area [-]
        dimensionedScalar fcl_;

        //- Tolerance criteria for iterative process to find Tcl
        scalar tolerance_;

        //- Maximum number of correctors for cloth temperature
        int maxClothIter_;

        //- Switch to use volume weighted velocity field for caluclation
        Switch meanVelocity_;


    // Private Member Functions

        // Calculate the magnitude of the velocity [m/s]
        tmp<volScalarField> magU() const;

        // Calculate the radiation temperature in the domain using a simple
        // approach [K]
        dimensionedScalar Trad() const;

        // Calculate the saturation pressure based on 7730:2006
        // Possible options: adding different calculation methods such as
        // the formulation based on Magnus or others [Pa]
        tmp<volScalarField> pSat() const;

        // Calculate and return the surface temperature of the cloth [K]
        // and the heat transfer coefficient hc [W/m^2/K]
        tmp<volScalarField> Tcloth
        (
            volScalarField& hc,
            const dimensionedScalar& metabolicRateSI,
            const dimensionedScalar& extWorkSI,
            const volScalarField& TdegC,
            const dimensionedScalar& Trad
        );

        // Return true if the cloth temperature iteration has converged
        bool converged(const volScalarField&) const;


public:

    //- Runtime type information
    TypeName("comfort");


    // Constructors

        //- Construct from Time and dictionary
        comfort
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~comfort();


    // Member Functions

        //- Read the data needed for the comfort calculation
        virtual bool read(const dictionary&);

        //- Calculate the predicted mean vote (PMV)
        //  and predicted percentage dissatisfaction (PPD) fields
        virtual bool execute();

        //- Write the PPD and PMV fields
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
