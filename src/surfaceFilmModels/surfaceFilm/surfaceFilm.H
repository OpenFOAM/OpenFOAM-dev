/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::regionModels::surfaceFilm

Description
    Base class for surface film models

SourceFiles
    surfaceFilm.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceFilm_H
#define surfaceFilm_H

#include "singleLayerRegionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{

/*---------------------------------------------------------------------------*\
                   Class surfaceFilm Declaration
\*---------------------------------------------------------------------------*/

class surfaceFilm
:
    public singleLayerRegionModel
{
    // Private data

        //- Acceleration due to gravity [m/s^2]
        const dimensionedVector& g_;


protected:

    // Protected member functions

        //- Read control parameters from dictionary
        virtual bool read();


public:

    //- Runtime type information
    TypeName("surfaceFilm");


    // Constructors

        //- Construct from type name, mesh and gravity vector
        surfaceFilm
        (
            const word& modelType,
            const fvMesh& mesh,
            const dimensionedVector& g,
            const word& regionType
        );

        //- Disallow default bitwise copy construction
        surfaceFilm(const surfaceFilm&) = delete;


    //- Destructor
    virtual ~surfaceFilm();


    // Member Functions

        // Access

            //- Return the acceleration due to gravity
            inline const dimensionedVector& g() const;

            //- External hook to add sources to the film
            virtual void addSources
            (
                const label patchi,
                const label facei,
                const scalar massSource,
                const vector& momentumSource,
                const scalar pressureSource,
                const scalar energySource
            ) = 0;


        // Fields

            //- Return the film thickness [m]
            virtual const volScalarField& delta() const = 0;

            //- Return the film coverage, 1 = covered, 0 = uncovered / []
            virtual const volScalarField& coverage() const = 0;

            //- Return the film velocity [m/s]
            virtual const volVectorField& U() const = 0;

            //- Return the film surface velocity [m/s]
            virtual tmp<volVectorField::Internal> Us() const = 0;

            //- Return the film wall velocity [m/s]
            virtual tmp<volVectorField::Internal> Uw() const = 0;

            //- Return the film density [kg/m^3]
            virtual const volScalarField& rho() const = 0;

            //- Return the film surface tension [N/m]
            virtual tmp<volScalarField> sigma() const = 0;


        // Transfer fields - to the primary region

            //- Return mass transfer source - Eulerian phase only
            virtual tmp<volScalarField> primaryMassTrans() const = 0;

            //- Return the film mass available for transfer
            virtual const volScalarField& cloudMassTrans() const = 0;

            //- Return the parcel diameters originating from film
            virtual const volScalarField& cloudDiameterTrans() const = 0;


        // Evolution

            //- Main driver routing to evolve the region - calls other evolves
            virtual void evolve();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const surfaceFilm&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "surfaceFilmI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
