/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::clouds::tracer

Description
    Cloud that follows the carrier flow. The elements have no size or mass and
    they have no physical modelling associated with them.

SourceFiles
    tracer.C

\*---------------------------------------------------------------------------*/

#ifndef tracer_H
#define tracer_H

#include "carried.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace clouds
{

/*---------------------------------------------------------------------------*\
                            Class tracer Declaration
\*---------------------------------------------------------------------------*/

class tracer
:
    public cloud,
    public carried
{
protected:

    // Protected Member Functions

        //- Return the acceleration with which to do second-order tracking
        virtual tmp<LagrangianSubVectorField> dUdt
        (
            const LagrangianSubMesh&
        ) const;

        //- Do we need to re-calculate particles that are modified?
        virtual bool reCalculateModified();

        //- Update the cloud properties
        virtual void calculate
        (
            const LagrangianSubScalarField& deltaT,
            const bool final
        );

        //- Partition hook
        virtual void partition();


public:

    // Public Static Data

        //- Run-time type information
        TypeName("tracer");


    // Constructors

        //- Construct from a mesh and context
        tracer
        (
            LagrangianMesh& mesh,
            const contextType context,
            const dictionary& dict
        );


    //- Destructor
    virtual ~tracer();


    // Member Functions

        //- Solve the cloud's evolution over the current time-step
        virtual void solve(const bool initial, const bool final);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace clouds
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
