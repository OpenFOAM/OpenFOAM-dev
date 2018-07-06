/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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
    Foam::waveSuperposition

Description
    A wrapper around a list of wave models. Superimposes the modelled values
    of elevation and velocity.

SourceFiles
    waveSuperposition.C

\*---------------------------------------------------------------------------*/

#ifndef waveSuperposition_H
#define waveSuperposition_H

#include "waveModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class waveSuperposition Declaration
\*---------------------------------------------------------------------------*/

class waveSuperposition
{
    // Private data

        //- Reference to the database
        const objectRegistry& db_;

        //- The origin of the wave coordinate system
        const vector origin_;

        //- The mean flow direction
        const vector direction_;

        //- The mean flow speed
        const scalar speed_;

        //- Wave models to superimpose
        PtrList<waveModel> waveModels_;

        //- The angle relative to the mean velocity at which the waves propagate
        scalarList waveAngles_;

        //- Ramp for the mean flow speed
        const autoPtr<Function1<scalar>> ramp_;

        //- Scaling in the flow direction
        const autoPtr<Function1<scalar>> scale_;

        //- Scaling perpendicular to the flow direction
        const autoPtr<Function1<scalar>> crossScale_;

        //- Calculate wave properties using the height above the wave (true) or
        //  the height above the origin (false)?
        const Switch heightAboveWave_;


    // Private Member Functions

        //- Get the transformation to actual coordinates
        void transformation
        (
            const vectorField& p,
            tensor& axes,
            scalar& u,
            vectorField& xyz
        ) const;

        //- Get the wave elevation relative to the mean at a given time, mean
        //  velocity and local coordinates. Local x is aligned with the mean
        //  velocity, and y is perpendicular to both x and gravity.
        tmp<scalarField> elevation
        (
            const scalar t,
            const vector2DField& xy
        ) const;

        //- Get the wave velocity at a given time, mean velocity and local
        //  coordinates. Local x is aligned with the mean velocity, z with
        //  negative gravity, and y is perpendicular to both.
        tmp<vectorField> velocity(const scalar t, const vectorField& xyz) const;

        //- Get the wave pressure at a given time, mean velocity and local
        //  coordinates. Local x is aligned with the mean velocity, z with
        //  negative gravity, and y is perpendicular to both.
        tmp<scalarField> pressure(const scalar t, const vectorField& xyz) const;

        //- Get the scaling factor, calculated from the optional scaling
        //  functions. X and y are the same as for the elevation method.
        tmp<scalarField> scale(const vector2DField& xy) const;


public:

    // Constructors

        //- Construct from a database
        waveSuperposition(const objectRegistry& db);

        //- Construct a copy
        waveSuperposition(const waveSuperposition& waves);

        //- Construct from a database and gravity
        waveSuperposition(const objectRegistry& db, const dictionary& dict);


    //- Destructor
    ~waveSuperposition();


    // Member Functions

        //- Get the height above the waves at a given time and global positions
        tmp<scalarField> height(const scalar t, const vectorField& p) const;

        //- Get the liquid velocity at a given time and global positions
        tmp<vectorField> ULiquid(const scalar t, const vectorField& p) const;

        //- Get the gas velocity at a given time and global positions
        tmp<vectorField> UGas(const scalar t, const vectorField& p) const;

        //- Get the liquid pressure at a given time and global positions
        tmp<scalarField> pLiquid(const scalar t, const vectorField& p) const;

        //- Get the gas pressure at a given time and global positions
        tmp<scalarField> pGas(const scalar t, const vectorField& p) const;

        //- Get the mean flow velocity
        inline vector UMean(const scalar t) const
        {
            return (ramp_.valid() ? ramp_->value(t) : 1)*direction_*speed_;
        }

        //- Write
        void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
