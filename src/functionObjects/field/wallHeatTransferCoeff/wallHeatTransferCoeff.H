/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2021 OpenFOAM Foundation
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
    Foam::functionObjects::wallHeatTransferCoeff

Description
    Calculates and writes the estimated heat transfer coefficient at wall
    patches as the volScalarField field.

    All wall patches are included by default; to restrict the calculation to
    certain patches, use the optional 'patches' entry.

    The models are selected run time by model entry. For detailed description
    look at the header file for specific model under
    wallHeatTransferCoeffModels.

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
        Pr          0.707;
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
        Pr          0.707;
        Prt         0.9;
        Lchar       0.001;
    }
    \endverbatim

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

Note
    Writing field 'wallHeatTransferCoeff' is done by default, but it can be
    overridden by defining an empty \c objects list. For details see
    writeLocalObjects.

See also
    Foam::functionObject
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::logFiles
    Foam::functionObjects::writeLocalObjects
    Foam::functionObjects::timeControl

SourceFiles
    wallHeatTransferCoeff.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_wallHeatTransferCoeff_H
#define functionObjects_wallHeatTransferCoeff_H

#include "fvMeshFunctionObject.H"
#include "logFiles.H"
#include "writeLocalObjects.H"
#include "wallHeatTransferCoeffModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                    Class wallHeatTransferCoeff Declaration
\*---------------------------------------------------------------------------*/

class wallHeatTransferCoeff
:
    public fvMeshFunctionObject,
    public logFiles,
    public writeLocalObjects
{

private:

    // Protected data

        // Pointer to the model
        autoPtr<wallHeatTransferCoeffModel> coeffModel_;

        // Return wall shear stress
        tmp<volSymmTensorField> devSigma();


    // Private data

        //- Density [kg/m^3]
        dimensionedScalar rho_;

        //- Specific capacity [J/K/kg]
        dimensionedScalar Cp_;


protected:

    // Protected data

        //- Reference to Time
        const Time& runTime_;

        //- Optional list of patches to process
        labelHashSet patchSet_;

        //- File header information
        virtual void writeFileHeader(const label i);


public:

    //- Runtime type information
    TypeName("wallHeatTransferCoeff");


    // Constructors

        //- Construct from name, mesh and dict
        wallHeatTransferCoeff
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~wallHeatTransferCoeff();


    // Member Functions

        //- Read the wallHeatTransferCoeffs data
        virtual bool read(const dictionary&);

        //- Calculate the wall heat transfer coefficient
        virtual bool execute();

        //- Write the wall heat transfer coefficient
        virtual bool write();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const wallHeatTransferCoeff&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
