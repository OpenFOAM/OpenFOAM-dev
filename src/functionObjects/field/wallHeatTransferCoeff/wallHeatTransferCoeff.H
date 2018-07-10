/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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
    Calculates and write the estimated incompressible flow heat transfer
    coefficient at wall patches as the volScalarField field
    'wallHeatTransferCoeff'.

    All wall patches are included by default; to restrict the calculation to
    certain patches, use the optional 'patches' entry.

    Example of function object specification:
    \verbatim
    wallHeatTransferCoeff1
    {
        type        wallHeatTransferCoeff;
        libs        ("libfieldFunctionObjects.so");
        ...
        region      fluid;
        patches     (".*Wall");
        rho         1.225;
        Cp          1005;
        Prl         0.707;
        Prt         0.9;
    }
    \endverbatim

Usage
    \table
        Property | Description                   | Required | Default value
        type     | Type name: wallHeatTransferCoeff | yes   |
        patches  | List of patches to process    | no       | all wall patches
        region   | Region to be evaluated        | no       | default region
        rho      | Fluid density                 | yes      |
        Cp       | Fluid heat capacity           | yes      |
        Prl      | Fluid laminar Prandtl number  | yes      |
        Prt      | Fluid turbulent Prandtl number| yes      |
    \endtable

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
#include "HashSet.H"
#include "volFieldsFwd.H"

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

protected:

    // Protected data

        //- Optional list of patches to process
        labelHashSet patchSet_;

        //- Fluid density
        scalar rho_;

        //- Fluid heat capacity
        scalar Cp_;

        //- Fluid laminar Prandtl number
        scalar Prl_;

        //- Fluid turbulent Prandtl number
        scalar Prt_;


    // Protected Member Functions

        //- File header information
        virtual void writeFileHeader(const label i);

        //- Calculate the heat transfer coefficient
        void calcHeatTransferCoeff
        (
            const volScalarField& nu,
            const volScalarField& nut,
            volScalarField& wallHeatTransferCoeff
        );


private:

    // Private member functions

        //- Disallow default bitwise copy construct
        wallHeatTransferCoeff(const wallHeatTransferCoeff&);

        //- Disallow default bitwise assignment
        void operator=(const wallHeatTransferCoeff&);


public:

    //- Runtime type information
    TypeName("wallHeatTransferCoeff");


    // Constructors

        //- Construct from Time and dictionary
        wallHeatTransferCoeff
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );


    //- Destructor
    virtual ~wallHeatTransferCoeff();


    // Member Functions

        //- Read the wallHeatTransferCoeff data
        virtual bool read(const dictionary&);

        //- Calculate the wall heat transfer coefficient
        virtual bool execute();

        //- Write the wall heat transfer coefficient
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
