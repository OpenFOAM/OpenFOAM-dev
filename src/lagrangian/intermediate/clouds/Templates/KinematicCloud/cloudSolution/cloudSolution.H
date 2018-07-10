/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::cloudSolution

Description
    Stores all relevant solution info for cloud

SourceFiles
    cloudSolutionI.H
    cloudSolution.C

\*---------------------------------------------------------------------------*/

#ifndef cloudSolution_H
#define cloudSolution_H

#include "fvMesh.H"
#include "Switch.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class cloudSolution Declaration
\*---------------------------------------------------------------------------*/

class cloudSolution
{
    // Private Data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Dictionary used during construction
        dictionary dict_;

        //- Cloud active flag
        const Switch active_;

        //- Transient flag
        Switch transient_;

        //- Calculation frequency - carrier steps per cloud step
        //  NOTE: Steady operation only
        label calcFrequency_;

        //- Maximum particle Courant number
        //  Max fraction of current cell that can be traversed in a single
        //  step
        scalar maxCo_;

        //- Current cloud iteration
        label iter_;

        //- Particle track time
        scalar trackTime_;


        // Run-time options

            //- Flag to indicate whether parcels are coupled to the carrier
            //  phase, i.e. whether or not to generate source terms for
            //  carrier phase
            Switch coupled_;

            //- Flag to correct cell values with latest transfer information
            //  during the lagrangian timestep
            Switch cellValueSourceCorrection_;

            //- Maximum particle track time [s]
            scalar maxTrackTime_;

            //- Flag to indicate whether coupling source terms should be
            //  reset on start-up/first read
            Switch resetSourcesOnStartup_;

            //- List schemes, e.g. U semiImplicit 1
            List<Tuple2<word, Tuple2<bool, scalar>>> schemes_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const cloudSolution&);


public:

    // Constructors

        //- Construct null from mesh reference
        cloudSolution(const fvMesh& mesh);

        //- Construct from mesh and dictionary
        cloudSolution(const fvMesh& mesh, const dictionary& dict);

        //- Construct copy
        cloudSolution(const cloudSolution& cs);


    //- Destructor
    virtual ~cloudSolution();


    // Member functions

        //- Read properties from dictionary
        void read();


        // Access

            //- Return relaxation coefficient for field
            scalar relaxCoeff(const word& fieldName) const;

            //- Return semi-implicit flag coefficient for field
            bool semiImplicit(const word& fieldName) const;

            //- Return reference to the mesh
            inline const fvMesh& mesh() const;

            //- Return const access to the dictionary
            inline const dictionary& dict() const;

            //- Return the active flag
            inline const Switch active() const;

            //- Return const access to the transient flag
            inline const Switch transient() const;

            //- Return const access to the steady flag
            inline const Switch steadyState() const;

            //- Return const access to the calculation frequency
            inline label calcFrequency() const;

            //- Return const access to the max particle Courant number
            inline scalar maxCo() const;

            //- Return const access to the current cloud iteration
            inline label iter() const;

            //- Increment and return iter counter
            inline label nextIter();

            //- Return the particle track time
            inline scalar trackTime() const;

            //- Return const access to the coupled flag
            inline const Switch coupled() const;

            //- Return non-const access to the coupled flag
            inline Switch& coupled();

            //- Return const access to the cell value correction flag
            inline const Switch cellValueSourceCorrection() const;

            //- Return const access to the particle track time
            inline scalar maxTrackTime() const;

            //- Return const access to the reset sources flag
            inline const Switch resetSourcesOnStartup() const;

            //- Source terms dictionary
            inline const dictionary& sourceTermDict() const;

            //- Interpolation schemes dictionary
            inline const dictionary& interpolationSchemes() const;

            //- Integration schemes dictionary
            inline const dictionary& integrationSchemes() const;


    // Helper functions

        //- Returns true if performing a cloud iteration this calc step
        bool solveThisStep() const;

        //- Returns true if possible to evolve the cloud and sets timestep
        //  parameters
        bool canEvolve();

        //- Returns true if writing this step
        bool output() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cloudSolutionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
