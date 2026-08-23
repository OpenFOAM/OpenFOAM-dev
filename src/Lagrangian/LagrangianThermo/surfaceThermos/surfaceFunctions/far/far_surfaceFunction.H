/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::surfaceFunctions::far

Description
    Surface function which uses only the far-field transport property value.
    For efficiency, this function should avoid evaluating the surface
    property entirely.

SourceFiles
    far_surfaceFunction.C

\*---------------------------------------------------------------------------*/

#ifndef far_surfaceFunction_H
#define far_surfaceFunction_H

#include "surfaceFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceFunctions
{

/*---------------------------------------------------------------------------*\
                          Class far Declaration
\*---------------------------------------------------------------------------*/

class far
:
    public surfaceFunction
{
    // Private Member Functions

        //- Return coefficients from which to construct a near-surface thermal
        //  conductivity
        virtual void thermalConductivityCoeffs
        (
            const LagrangianSubMesh& subMesh,
            autoPtr<dimensionedScalar>& fValuePtr,
            tmp<LagrangianSubScalarField>& tfField
        ) const;

        //- Return coefficients from which to construct a near-surface dynamic
        //  diffusivity or viscosity
        virtual void dynamicDiffusivityCoeffs
        (
            const LagrangianSubMesh& subMesh,
            autoPtr<dimensionedScalar>& fValuePtr,
            tmp<LagrangianSubScalarField>& tfField
        ) const;

        //- Return coefficients from which to construct a near-surface kinematic
        //  diffusivity or viscosity
        virtual void kinematicDiffusivityCoeffs
        (
            const LagrangianSubMesh& subMesh,
            autoPtr<dimensionedScalar>& fValuePtr,
            tmp<LagrangianSubScalarField>& tfField
        ) const;


public:

    //- Runtime type information
    TypeName("far");


    // Constructors

        //- Construct for a given mesh
        far(const LagrangianMesh& mesh);


    //- Destructor
    virtual ~far();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
