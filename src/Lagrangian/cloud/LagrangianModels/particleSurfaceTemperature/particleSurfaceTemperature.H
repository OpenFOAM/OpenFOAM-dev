/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025-2026 OpenFOAM Foundation
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
    Foam::Lagrangian::particleSurfaceTemperature

Description
    Model which sets the surface temperature to be the same as the internal
    temperature of the particles. This implies that the heat transfer processes
    within the particles are fast in comparison with those outside. This is
    typically a reasonable assumption for solid particles or liquid droplets
    surrounded by a relatively insulating gaseous carrier.

Usage
    Example specification:
    \verbatim
    <LagrangianModelName>
    {
        type        particleSurfaceTemperature;
    }
    \endverbatim

SourceFiles
    particleSurfaceTemperature.C

\*---------------------------------------------------------------------------*/

#ifndef particleSurfaceTemperature_H
#define particleSurfaceTemperature_H

#include "cloudLagrangianModel.H"
#include "LagrangianSource.H"
#include "sharedRegIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Lagrangian
{

/*---------------------------------------------------------------------------*\
                    Class particleSurfaceTemperature Declaration
\*---------------------------------------------------------------------------*/

class particleSurfaceTemperature
:
    public LagrangianModel,
    private cloudLagrangianModel
{
public:

    //- Runtime type information
    TypeName("particleSurfaceTemperature");


    // Constructors

        //- Construct from components
        particleSurfaceTemperature
        (
            const word& name,
            const LagrangianMesh& mesh,
            const dictionary& modelDict,
            const dictionary& stateDict
        );


    // Member Functions

        //- Return false as this model does not create sources
        virtual bool addsSupToField(const word&, const word&) const;

        //- Update the surface temperature
        virtual void calculate
        (
            const LagrangianSubScalarField& deltaT,
            const bool final
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Lagrangian
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
