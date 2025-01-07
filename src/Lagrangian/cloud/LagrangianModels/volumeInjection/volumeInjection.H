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
    Foam::volumeInjection

Description
    Volume injection model. This injects particles instantaneously within the
    volume of a set of cells. Particles are positioned randomly.

Usage
    \table
        Property      | Description                    | Required? | Default
        select        | Type of cell sub-set to \
                        inject into, or 'all'          | maybe     | <inferred>
        cellSet       | The name of the cell set       | if select is cellSet |
        cellZone      | The name of the cell zone      | if select is cellZone |
        points        | Points inside the cells        | if select is points |
        number        | The number to inject           | if numberDensity is \
                                                         not specified |
        numberDensity | The number to inject per unit \
                        volume                         | if number is not \
                                                         specified |
        time          | The time at which to inject    | no        | 0
    \endtable

    Example specification:
    \verbatim
    <LagrangianModelName>
    {
        type        volumeInjection;
        cellZone    bed;
        numberDensity 2 [mm^-3];
        time        0.01;
    }
    \endverbatim

See also
    Foam::polyCellSet

SourceFiles
    volumeInjection.C

\*---------------------------------------------------------------------------*/

#ifndef volumeInjection_H
#define volumeInjection_H

#include "LagrangianInjection.H"
#include "polyCellSet.H"
#include "restartableRandomGenerator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Lagrangian
{

/*---------------------------------------------------------------------------*\
                          Class volumeInjection Declaration
\*---------------------------------------------------------------------------*/

class volumeInjection
:
    public LagrangianInjection
{
private:

    // Private Data

        //- The set of cells into which to inject
        mutable polyCellSet set_;

        //- Do we have the number (true) or the number density (false)
        bool haveNumber_;

        //- The number to inject or the number density with which to inject
        scalar numberOrNumberDensity_;

        //- The time at which to inject
        scalar time_;

        //- A global random number generator
        mutable restartableRandomGenerator globalRndGen_;

        //- A local random number generator
        mutable restartableRandomGenerator localRndGen_;

        //- The time index
        mutable label timeIndex_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs(const dictionary& modelDict);


public:

    //- Runtime type information
    TypeName("volumeInjection");


    // Constructors

        //- Construct from components
        volumeInjection
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
