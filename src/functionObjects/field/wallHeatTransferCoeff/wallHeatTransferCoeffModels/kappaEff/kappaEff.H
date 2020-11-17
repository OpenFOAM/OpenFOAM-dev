/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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
    Foam::wallHeatTransferCoeffModels::kappaEff

Description
    Calculates the estimated flow heat transfer coefficient at wall patches
    as the volScalarField field 'kappaEff' using one of equeations bellow.

    kappaEff model, given by:

    \f[
        htc = \rho*C_p*({\frac{\nu}{Prl} + \frac{\nu_t}{Prt}})
    \f]

    kappaEff model with characteristic length, given by:

    \f[
        htc =
            \rho*C_p*({\frac{\nu}{Prl} + \frac{\nu_t}{Prt}})\frac{1}{L_{char}};
    \f]

    where
    \vartable
        rho      | Density [kg/m^3]
        Cp       | Specific heat capacity [m^2/K/s^2)]
        Prl      | Fluid laminar Prandtl number []
        Prt      | Fluid turbulent Prandtl number []
        Lchar    | Characteristic length [m]
    \endvartable

    Example of function object specification:
    \verbatim
    kappaEff1
    {
        type        wallHeatTransferCoeff;
        libs        ("libfieldFunctionObjects.so");
        model       kappaEff;
        ...
        region      fluid;
        patches     (".*Wall");
        rho         1.225;
        Cp          1005;
        Prl         0.707;
        Prt         0.9;
    }
    \endverbatim

    \verbatim
    kappaEff2
    {
        type        wallHeatTransferCoeff;
        libs        ("libfieldFunctionObjects.so");
        model       kappaEff;
        ...
        region      fluid;
        patches     (".*Wall");
        rho         1.225;
        Cp          1005;
        Prl         0.707;
        Prt         0.9;
        Lchar       0.001;
    }
    \endverbatim

Usage
    \table
        Property | Description                   | Required | Default value
        type     | Type name: wallHeatTransferCoeff | yes   |
        model    | Type name: kappaEff           | no       | kappaEff
        patches  | List of patches to process    | no       | all wall patches
        region   | Region to be evaluated        | no       | default region
        rho      | Fluid density                 | yes      |
        Cp       | Fluid heat capacity           | yes      |
        Prl      | Fluid laminar Prandtl number  | yes      |
        Prt      | Fluid turbulent Prandtl number| yes      |
        Lchar    | Characteristic length         | no       | no
    \endtable

Note
    Cp and rho are required only for incompressible flow calclulations.

SourceFiles
    kappaEff.C

\*---------------------------------------------------------------------------*/

#ifndef wallHeatTransferCoeffModels_kappaEff_H
#define wallHeatTransferCoeffModels_kappaEff_H

#include "wallHeatTransferCoeffModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallHeatTransferCoeffModels
{

/*---------------------------------------------------------------------------*\
                          Class kappaEff Declaration
\*---------------------------------------------------------------------------*/

class kappaEff
:
    public wallHeatTransferCoeffModel
{
    // Private data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Fluid laminar Prandtl number
        dimensionedScalar Prl_;

        //- Fluid turbulent Prandtl number
        dimensionedScalar Prt_;

        //- Characteristic length
        dimensionedScalar Lchar_;

        //- Is characteristic length used?
        bool isCharLength_;


public:

    //- Runtime type information
    TypeName("kappaEff");


    // Constructors

        //- Construct from name, mesh and dict
        kappaEff
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        kappaEff(const kappaEff&) = delete;


    //- Destructor
    virtual ~kappaEff();


    // Member Functions

        //- Read the kappaEff data
        virtual bool read(const dictionary&);

        //- Calculate the wall heat transfer coefficient
        virtual tmp<volScalarField> htcByRhoCp
        (
            const momentumTransportModel& mmtm,
            const labelHashSet& patches
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const kappaEff&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace wallHeatTransferCoeffModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
