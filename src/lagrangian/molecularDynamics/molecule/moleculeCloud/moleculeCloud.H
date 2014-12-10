/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::moleculeCloud

Description

SourceFiles
    moleculeCloudI.H
    moleculeCloud.C


\*---------------------------------------------------------------------------*/

#ifndef moleculeCloud_H
#define moleculeCloud_H

#include "Cloud.H"
#include "molecule.H"
#include "IOdictionary.H"
#include "potential.H"
#include "InteractionLists.H"
#include "labelVector.H"
#include "Random.H"
#include "fileName.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class moleculeCloud Declaration
\*---------------------------------------------------------------------------*/

class moleculeCloud
:
    public Cloud<molecule>
{

private:

    // Private data

        const polyMesh& mesh_;

        const potential& pot_;

        List<DynamicList<molecule*> > cellOccupancy_;

        InteractionLists<molecule> il_;

        List<molecule::constantProperties> constPropList_;

        Random rndGen_;


    // Private Member Functions

        void buildConstProps();

        void setSiteSizesAndPositions();

        //- Determine which molecules are in which cells
        void buildCellOccupancy();

        void calculatePairForce();

        inline void evaluatePair
        (
            molecule& molI,
            molecule& molJ
        );

        inline bool evaluatePotentialLimit
        (
            molecule& molI,
            molecule& molJ
        ) const;

        void calculateTetherForce();

        void calculateExternalForce();

        void removeHighEnergyOverlaps();

        void initialiseMolecules
        (
            const IOdictionary& mdInitialiseDict
        );

        void createMolecule
        (
            const point& position,
            label cell,
            label tetFace,
            label tetPt,
            label id,
            bool tethered,
            scalar temperature,
            const vector& bulkVelocity
        );

        label nSites() const;

        inline vector equipartitionLinearVelocity
        (
            scalar temperature,
            scalar mass
        );

        inline vector equipartitionAngularMomentum
        (
            scalar temperature,
            const molecule::constantProperties& cP
        );

        //- Disallow default bitwise copy construct
        moleculeCloud(const moleculeCloud&);

        //- Disallow default bitwise assignment
        void operator=(const moleculeCloud&);


public:

    // Constructors

        //- Construct given mesh and potential references
        moleculeCloud
        (
            const polyMesh& mesh,
            const potential& pot,
            bool readFields = true
        );

        //- Construct given mesh, potential and mdInitialiseDict
        moleculeCloud
        (
            const polyMesh& mesh,
            const potential& pot,
            const IOdictionary& mdInitialiseDict,
            bool readFields = true
        );


    // Member Functions

        //- Evolve the molecules (move, calculate forces, control state etc)
        void evolve();

        void calculateForce();

        void applyConstraintsAndThermostats
        (
            const scalar targetTemperature,
            const scalar measuredTemperature
        );


        // Access

            inline const polyMesh& mesh() const;

            inline const potential& pot() const;

            inline const List<DynamicList<molecule*> >& cellOccupancy() const;

            inline const InteractionLists<molecule>& il() const;

            inline const List<molecule::constantProperties> constProps() const;

            inline const molecule::constantProperties&
                constProps(label id) const;

            inline Random& rndGen();


    // Member Operators

        //- Write molecule sites in XYZ format
        void writeXYZ(const fileName& fName) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "moleculeCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
