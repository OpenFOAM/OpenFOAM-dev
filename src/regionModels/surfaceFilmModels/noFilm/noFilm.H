/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::noFilm

Description
    Dummy surface film model for 'none'

SourceFiles
    noFilm.C

\*---------------------------------------------------------------------------*/

#ifndef noFilm_H
#define noFilm_H

#include "surfaceFilmModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                          Class noFilm Declaration
\*---------------------------------------------------------------------------*/

class noFilm
:
    public surfaceFilmModel
{
private:

    // Private member functions

        //- Disallow default bitwise copy construct
        noFilm(const noFilm&);

        //- Disallow default bitwise assignment
        void operator=(const noFilm&);


protected:

    // Protected member functions

        //- Read control parameters from dictionary
        virtual bool read();


public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from components
        noFilm
        (
            const word& modelType,
            const fvMesh& mesh,
            const dimensionedVector& g,
            const word& regionType
        );


    //- Destructor
    virtual ~noFilm();


    // Member Functions

        // Access

            //- External hook to add sources to the film
            virtual void addSources
            (
                const label patchI,
                const label faceI,
                const scalar massSource,
                const vector& momentumSource,
                const scalar pressureSource,
                const scalar energySource
            );


        // Fields

            //- Return the film thickness [m]
            virtual const volScalarField& delta() const;

            //- Return the film coverage, 1 = covered, 0 = uncovered / []
            virtual const volScalarField& alpha() const;

            //- Return the film velocity [m/s]
            virtual const volVectorField& U() const;

            //- Return the film density [kg/m3]
            virtual const volScalarField& rho() const;

            //- Return the film surface velocity [m/s]
            virtual const volVectorField& Us() const;

            //- Return the film wall velocity [m/s]
            virtual const volVectorField& Uw() const;

            //- Return the film mean temperature [K]
            virtual const volScalarField& T() const;

            //- Return the film surface temperature [K]
            virtual const volScalarField& Ts() const;

            //- Return the film wall temperature [K]
            virtual const volScalarField& Tw() const;

            //- Return the film specific heat capacity [J/kg/K]
            virtual const volScalarField& Cp() const;

            //- Return the film thermal conductivity [W/m/K]
            virtual const volScalarField& kappa() const;

            //- Return const access to the surface tension / [m/s2]
            inline const volScalarField& sigma() const;


            // Transfer fields - to the primary region

                //- Return mass transfer source - Eulerian phase only
                virtual tmp<volScalarField> primaryMassTrans() const;

                //- Return the film mass available for transfer
                virtual const volScalarField& cloudMassTrans() const;

                //- Return the parcel diameters originating from film
                virtual const volScalarField& cloudDiameterTrans() const;


        // Source fields

            // Mapped into primary region

                //- Return total mass source - Eulerian phase only
                virtual tmp<DimensionedField<scalar, volMesh> > Srho() const;

                //- Return mass source for specie i - Eulerian phase only
                virtual tmp<DimensionedField<scalar, volMesh> > Srho
                (
                    const label i
                ) const;

                //- Return enthalpy source - Eulerian phase only
                virtual tmp<DimensionedField<scalar, volMesh> > Sh() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
