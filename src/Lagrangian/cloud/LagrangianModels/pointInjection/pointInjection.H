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
    Foam::pointInjection

Description
    Point injection model. This injects particles continuously at a point with
    a given number rate. The point and the number rate are both Function1-s and
    can both vary with time.

    Note that this model only controls the number and position of injected
    Lagrangian particles. All physical properties are specified by
    corresponding source conditions. So the velocity/direction/angle/etc..., is
    controlled by the velocity source condition, the size distribution by the
    diameter source condition, and the flow rate by the number source
    condition.

Usage
    \table
        Property    | Description                     | Required? | Default
        point       | The position at which to inject | yes       |
        numberRate  | The number of particles to inject per unit time | yes |
    \endtable

    Example specification:
    \verbatim
    <LagrangianModelName>
    {
        type        pointInjection;
        point       (1 2 3) [cm];
        numberRate  10000;
    }
    \endverbatim

See also
    Foam::Function1s

SourceFiles
    pointInjection.C

\*---------------------------------------------------------------------------*/

#ifndef pointInjection_H
#define pointInjection_H

#include "LagrangianInjection.H"
#include "Function1.H"
#include "restartableRandomGenerator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Lagrangian
{

/*---------------------------------------------------------------------------*\
                          Class pointInjection Declaration
\*---------------------------------------------------------------------------*/

class pointInjection
:
    public LagrangianInjection
{
private:

    // Private Data

        //- The point at which to inject
        autoPtr<Function1<point>> point_;

        //- The rate at which to inject
        autoPtr<Function1<scalar>> numberRate_;

        //- The number deferred to the next injection step
        mutable scalar numberDeferred_;

        //- A random number generator
        mutable restartableRandomGenerator rndGen_;

        //- The time index
        mutable label timeIndex_;

        //- Where is the injection?
        enum class injectionLocation
        {
            unset,
            fixedPointOnThisProcessor,
            fixedPointOnAnotherProcessor,
            multiplePoints
        } injectionLocation_;

        //- The coordinates of the point
        barycentric coordinates_;

        //- The cell containing the point
        label celli_;

        //- The face adjacent to the tet containing the point
        label facei_;

        //- The face-tri adjacent to the tet containing the point
        label faceTrii_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs(const dictionary& modelDict);


public:

    //- Runtime type information
    TypeName("pointInjection");


    // Constructors

        //- Construct from components
        pointInjection
        (
            const word& name,
            const LagrangianMesh& mesh,
            const dictionary& modelDict,
            const dictionary& stateDict
        );


    // Member Functions

        //- Correct the LagrangianModel
        virtual void correct();

        //- Create new elements in the Lagrangian mesh
        LagrangianSubMesh modify
        (
            LagrangianMesh& mesh,
            const LagrangianSubMesh&
        ) const;


        // Mesh changes

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&);


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& modelDict);

            //- Write state
            virtual void writeState(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Lagrangian
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
