/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2020 OpenFOAM Foundation
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
    A wrapper around a list of wave models. Superimposes the modelled values of
    elevation and velocity. The New method looks up or or constructs an
    instance of this class on demand and returns a reference. Properties are
    read from a dictionary in constant.

Usage
    \table
        Property   | Description                             | Req'd? | Default
        origin     | origin of the wave coordinate system    | yes    |
        direction  | direction of the wave coordinate system | yes    |
        waves      | list of wave models to superimpose      | yes    |
        UMean      | velocity of the mean flow               | yes    |
        scale      | scale factor in the direction           | no     | None
        crossScale | scale factor perpendicular to the direction | no | None
        heightAboveWave | use the height above the wave as the vertical \\
        coordinate | no | false
    \endtable

    Example specification:
    \verbatim
    origin      (0 25 0);
    direction   (1 0 0);
    waves
    (
        Airy
        {
            length      40;
            amplitude   0.5;
            phase       0;
            angle       0;
        }
        Airy
        {
            length      20;
            amplitude   0.25;
            phase       1.5708;
            angle       0;
        }
    );
    UMean       (2 0 0);
    scale       table ((100 1) (200 0));
    crossScale  constant 1;
    heightAboveWave no;
    \endverbatim

SourceFiles
    waveSuperposition.C
    waveSuperpositionNew.C

\*---------------------------------------------------------------------------*/

#ifndef waveSuperposition_H
#define waveSuperposition_H

#include "waveModel.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class waveSuperposition Declaration
\*---------------------------------------------------------------------------*/

class waveSuperposition
:
    public IOdictionary
{
protected:

    // Protected Data

        //- The origin of the wave coordinate system
        const vector origin_;

        //- The direction of the wave coordinate system
        const vector direction_;

        //- Wave models to superimpose
        PtrList<waveModel> waveModels_;

        //- The angle relative to the direction at which the waves propagate
        scalarList waveAngles_;

        //- Mean velocity
        const autoPtr<Function1<vector>> UMean_;

        //- Scaling in the local x-direction
        const autoPtr<Function1<scalar>> scale_;

        //- Scaling perpendicular to the local x-direction
        const autoPtr<Function1<scalar>> crossScale_;

        //- Calculate wave properties using the height above the wave (true) or
        //  the height above the origin (false)?
        const Switch heightAboveWave_;


    // Protected Member Functions

        //- Get the transformation to actual coordinates
        void transformation
        (
            const scalar t,
            const vectorField& p,
            tensor& axes,
            vectorField& xyz
        ) const;

        //- Get the wave elevation relative to the mean at a given time and
        //  local coordinates. Local x is aligned with the direction, and y is
        //  perpendicular to both x and gravity.
        tmp<scalarField> elevation
        (
            const scalar t,
            const vector2DField& xy
        ) const;

        //- Get the wave velocity at a given time and local coordinates. Local
        //  x is aligned with the direction, z with negative gravity, and y is
        //  perpendicular to both.
        tmp<vectorField> velocity(const scalar t, const vectorField& xyz) const;

        //- Get the wave pressure at a given time and local coordinates. Local
        //  x is aligned with the direction, z with negative gravity, and y is
        //  perpendicular to both.
        tmp<scalarField> pressure(const scalar t, const vectorField& xyz) const;

        //- Get the scaling factor, calculated from the optional scaling
        //  functions. X and y are the same as for the elevation method.
        tmp<scalarField> scale(const vector2DField& xy) const;


public:

    //- Runtime type information
    TypeName("wave");


    // Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        waveSuperposition,
        objectRegistry,
        (const objectRegistry& db),
        (db)
    );


    // Static Data

        //- The name of the dictionary
        static const word dictName;


    // Static Member Functions

        //- Return a reference to the wave model on the given database,
        //  constructing if it doesn't exist
        static const waveSuperposition& New(const objectRegistry& db);


    // Constructors

        //- Construct from a database
        waveSuperposition(const objectRegistry& db);


    //- Destructor
    ~waveSuperposition();


    // Member Functions

        //- Get the height above the waves at a given time and global positions
        virtual tmp<scalarField> height
        (
            const scalar t,
            const vectorField& p
        ) const;

        //- Get the liquid velocity at a given time and global positions
        virtual tmp<vectorField> ULiquid
        (
            const scalar t,
            const vectorField& p
        ) const;

        //- Get the gas velocity at a given time and global positions
        virtual tmp<vectorField> UGas
        (
            const scalar t,
            const vectorField& p
        ) const;

        //- Get the liquid pressure at a given time and global positions
        virtual tmp<scalarField> pLiquid
        (
            const scalar t,
            const vectorField& p
        ) const;

        //- Get the gas pressure at a given time and global positions
        virtual tmp<scalarField> pGas
        (
            const scalar t,
            const vectorField& p
        ) const;

        //- Inherit write from regIOobject
        using regIOobject::write;

        //- Write
        void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
