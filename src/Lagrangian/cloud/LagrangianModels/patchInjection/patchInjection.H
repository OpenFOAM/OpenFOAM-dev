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
    Foam::patchInjection

Description
    Patch injection model. This injects particles continuously at a patch with
    a given number rate. The number rate is a Function1 and can vary with time.

    Note that this model only controls the number and position of injected
    Lagrangian particles. All physical properties are specified by
    corresponding source conditions. So the velocity/direction/angle/etc..., is
    controlled by the velocity source condition, the size distribution by the
    diameter source condition, and the flow rate by the number source
    condition.

Usage
    \table
        Property    | Description                     | Required? | Default
        patch       | Name of the patch               | yes       |
        numberRate  | The number of particles to inject per unit time | yes |
    \endtable

    Example specification:
    \verbatim
    <LagrangianModelName>
    {
        type        patchInjection;
        patch       top;
        numberRate  100000;
    }
    \endverbatim

See also
    Foam::Function1s

SourceFiles
    patchInjection.C

\*---------------------------------------------------------------------------*/

#ifndef patchInjection_H
#define patchInjection_H

#include "LagrangianInjection.H"
#include "cloudLagrangianModel.H"
#include "Function1.H"
#include "restartableRandomGenerator.H"
#include "CompactListList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Lagrangian
{

/*---------------------------------------------------------------------------*\
                          Class patchInjection Declaration
\*---------------------------------------------------------------------------*/

class patchInjection
:
    public LagrangianInjection,
    private cloudLagrangianModel
{
private:

    // Private Data

        //- The name of the patch through which to inject
        word patchName_;

        //- The rate at which to inject
        autoPtr<Function1<scalar>> numberRate_;

        //- The number deferred to the next injection step
        mutable scalar numberDeferred_;

        //- A global random number generator
        mutable restartableRandomGenerator globalRndGen_;

        //- A local random number generator
        mutable restartableRandomGenerator localRndGen_;

        //- The time index
        mutable label timeIndex_;

        //- Cumulative sum of the total patch areas in each process
        scalarList procSumArea_;

        //- Cumulative sum of the face areas in the patch
        scalarList patchFaceSumArea_;

        //- Cumulative sum of the triangle areas in each patch face
        CompactListList<scalar> patchFaceTriSumArea_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs(const dictionary& modelDict);

        //- Initialise the area sums
        void calcSumAreas();


public:

    //- Runtime type information
    TypeName("patchInjection");


    // Constructors

        //- Construct from components
        patchInjection
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

            //- Write state
            virtual void writeProcessorState(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Lagrangian
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
