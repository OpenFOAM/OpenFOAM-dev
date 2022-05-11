/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2022 OpenFOAM Foundation
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
    Foam::wallHeatTransferCoeffModels::ReynoldsAnalogy

Description
    Calculates and writes the estimated flow heat transfer coefficient at wall
    patches as the volScalarField field 'wallHeatTransferCoeff' using Reynolds
    Analogy.

    Reynolds Analogy model, given by:

    \f[
        htc = 0.5 \rho C_p |U_{ref}| C_f
    \f]

    where
    \vartable
        rho       | Density [kg/m^3]
        Cp        | Specific heat capacity [m^2/K/s^2)]
        Uref      | Far field velocity magnitude [m/s]
        Cf        | Skin friction coefficient []
    \endvartable


    All wall patches are included by default; to restrict the calculation to
    certain patches, use the optional 'patches' entry.

    Example of function object specification:
    \verbatim
    ReynoldsAnalogy1
    {
        type       wallHeatTransferCoeff;
        libs       ("libfieldFunctionObjects.so");
        model      ReynoldsAnalogy;
        ...
        region     fluid;
        patches    (".*Wall");
        rho        1.225;
        Cp         1005;
        Uref       1.0;
    }
    \endverbatim

Usage
    \table
        Property | Description               | Required | Default value
        type     | type name: wallHeatTransferCoeff | yes   |
        model    | Type name: ReynoldsAnalogy    | no       | kappaEff
        patches  | List of patches to process    | no       | all wall patches
        region   | Region to be evaluated        | no       | default region
        rho      | Density                       | yes      |
        Cp       | Specific heat capacity        | yes      |
        Uref     | Reference velocity magnitude  | yes      | no
    \endtable

    Note:
        Cp and rho are required only for incompressible flow calclulations.

SourceFiles
    ReynoldsAnalogy.C

\*---------------------------------------------------------------------------*/

#ifndef wallHeatTransferCoeffModels_ReynoldsAnalogy_H
#define wallHeatTransferCoeffModels_ReynoldsAnalogy_H

#include "wallHeatTransferCoeffModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallHeatTransferCoeffModels
{

/*---------------------------------------------------------------------------*\
                       Class ReynoldsAnalogy Declaration
\*---------------------------------------------------------------------------*/

class ReynoldsAnalogy
:
    public wallHeatTransferCoeffModel
{
    // Private Data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Reference velocity magnitude [m/s]
        dimensionedScalar Uref_;


    // Private member functions

        //- Calculate the wall shear-stress
        tmp<volVectorField> calcShearStress
        (
            const volSymmTensorField& tau
        );


public:

    //- Runtime type information
    TypeName("ReynoldsAnalogy");


    // Constructors

        //- Construct from name, mesh and dict
        ReynoldsAnalogy
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        ReynoldsAnalogy
        (
            const ReynoldsAnalogy&
        ) = delete;


    //- Destructor
    virtual ~ReynoldsAnalogy();


    // Member Functions

        //- Read the ReynoldsAnalogy data
        virtual bool read(const dictionary&);

        //- Calculate the heat transfer coefficient
        virtual tmp<volScalarField> htcByRhoCp
        (
            const momentumTransportModel& mmtm,
            const labelHashSet& patches
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ReynoldsAnalogy&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace wallHeatTransferCoeffModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
