/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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
    Foam::solvers::buoyancy

Description
    Buoyancy related data for the Foam::solvers::isothermalFluid solver module
    when solving buoyant cases with \c p_rgh and is selected based on the
    presence of the \c p_rgh field file.

    Provides:
        g : gravitational acceleration
        hRef : optional reference height
        pRef : optional reference pressure
        ghRef : -mag(g)*hRef
        gh : (g & h) - ghRef
        ghf : (g & hf) - ghRef
        p_rgh : p - rho*gh - pRef

    The \c Foam::solvers::buoyancy::New function returns an \c autoPtr<buoyancy>
    pointer containing either a \c buoyancy class pointer if the \c p_rgh field
    file is available otherwise a \c nullptr which can be checked using
    the \c autoPtr::valid() member function.

SourceFiles
    buoyancy.C

\*---------------------------------------------------------------------------*/

#ifndef buoyancy_H
#define buoyancy_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class buoyancy Declaration
\*---------------------------------------------------------------------------*/

class buoyancy
{
    // Private member data

        const fvMesh& mesh;

        const Time& runTime;


public:

    // Thermophysical properties

        //- Gravitational acceleration
        uniformDimensionedVectorField g;

        //- Reference height
        uniformDimensionedScalarField hRef;

        //- Reference pressure
        uniformDimensionedScalarField pRef;

        //- -mag(g)*hRef
        dimensionedScalar ghRef;

        //- (g & h) - ghRef
        volScalarField gh;

        //- (g & hf) - ghRef
        surfaceScalarField ghf;

        //- Buoyant pressure p - rho*gh - pRef
        volScalarField p_rgh;


    // Constructors

        //- Construct from the region mesh
        buoyancy(const fvMesh& mesh);

        //- Disallow default bitwise copy construction
        buoyancy(const buoyancy&) = delete;


    // Selectors

        //- Select, construct and return the buoyancy
        static autoPtr<buoyancy> New(const fvMesh& mesh);


    //- Destructor
    virtual ~buoyancy();


    // Member Functions

        //- Update gh and ghf following mesh-motion
        void moveMesh();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const buoyancy&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
