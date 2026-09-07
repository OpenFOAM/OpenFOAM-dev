/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2026 OpenFOAM Foundation
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
    Foam::fv::limitGroupSurfaceAreaVolumeRatio

Description
    Limits each group surface-area-volume-ratio (kappa) field within a
    population balance that uses a fractal shape model such that the value is
    greater than that of spherical particles represented by the corresponding
    group, but less than that of a minimum representative diameter.

    The minimum representative diameter can be specified by the user, but if
    omitted will default to the representative diameter of the first group in
    the population balance.

Usage
    Example usage:
    \verbatim
    limitGroupSurfaceAreaVolumeRatio1
    {
        type            limitGroupSurfaceAreaVolumeRatio;
        libs            ("libmultiphaseEulerFvConstraints.so");

        populationBalance  aggregates;

        // Optional: Minimum diameter defaults to the representative diameter
        // of the first group in the population balance
        minDiameter        4e-10;
    }
    \endverbatim

SourceFiles
    limitGroupSurfaceAreaVolumeRatio.C

\*---------------------------------------------------------------------------*/

#ifndef limitGroupSurfaceAreaVolumeRatio_H
#define limitGroupSurfaceAreaVolumeRatio_H

#include "fvConstraint.H"
#include "populationBalanceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
              Class limitGroupSurfaceAreaVolumeRatio Declaration
\*---------------------------------------------------------------------------*/

class limitGroupSurfaceAreaVolumeRatio
:
    public fvConstraint
{
    // Private Data

        //- Population balance model
        const populationBalanceModel& popBal_;

        //- Minimum diameter used to define the upper kappa bound
        dimensionedScalar minDiameter_;

        //- List of constrained field names
        wordList constrainedFields_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs(const dictionary& dict);


public:

    //- Runtime type information
    TypeName("limitGroupSurfaceAreaVolumeRatio");


    // Constructors

        //- Construct from components
        limitGroupSurfaceAreaVolumeRatio
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        limitGroupSurfaceAreaVolumeRatio
        (
            const limitGroupSurfaceAreaVolumeRatio&
        ) = delete;


    //- Destructor
    virtual ~limitGroupSurfaceAreaVolumeRatio()
    {}


    // Member Functions

        //- Return the list of fields constrained by the fvConstraint
        virtual wordList constrainedFields() const;

        //- Apply constraint to the kappa field
        virtual bool constrain(VolField<scalar>& field) const;

        //- Update for mesh motion
        virtual bool movePoints();

        //- Update topology using the given map
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);

        //- Redistribute or update using the given distribution map
        virtual void distribute(const polyDistributionMap&);

        //- Read dictionary
        virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const limitGroupSurfaceAreaVolumeRatio&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
