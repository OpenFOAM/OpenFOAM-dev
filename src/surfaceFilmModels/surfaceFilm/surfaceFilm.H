/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
    Foam::surfaceFilm

Description
    Base class for surface film models

SourceFiles
    surfaceFilm.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceFilm_H
#define surfaceFilm_H

#include "fvMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class surfaceFilm Declaration
\*---------------------------------------------------------------------------*/

class surfaceFilm
:
    public IOdictionary
{
protected:

    // Protected data

        //- Reference to the primary mesh
        const fvMesh& primaryMesh_;

        //- Reference to the time
        const Time& time_;

        //- Information output
        Switch infoOutput_;

        //- Model type
        const word modelType_;

        //- Film name
        word regionName_;

        //- Film mesh
        fvMesh mesh_;

        //- Model coefficients dictionary
        dictionary coeffs_;

        //- Dictionary of output properties
        autoPtr<IOdictionary> outputPropertiesPtr_;


        // Addressing

            //- List of patch IDs on the primary film coupled to this film
            labelList primaryPatchIDs_;

            //- List of patch IDs internally coupled with the primary film
            labelList intCoupledPatchIDs_;

            //- List of patch IDs opposite to internally coupled patches
            labelList passivePatchIDs_;


        //- Patch normal vectors
        volVectorField nHat_;

        //- Cell cross-sectional area magnitude [m^2]
        volScalarField::Internal magSf_;

        //- Cell layer volume/area [m]
        volScalarField VbyA_;

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
            const word& filmType
        );

        //- Disallow default bitwise copy construction
        surfaceFilm(const surfaceFilm&) = delete;


    //- Destructor
    virtual ~surfaceFilm();


    // Member Functions

        // Access

            //- Return the reference to the time
            inline const Time& time() const;

            //- Return the reference to the primary mesh
            inline const fvMesh& primaryMesh() const;

            //- Return the film mesh
            inline const fvMesh& mesh() const;

            //- Return the information flag
            inline const Switch& infoOutput() const;

            //- Return the model coefficients dictionary
            inline const dictionary& coeffs() const;

            //- Return the solution dictionary
            inline const dictionary& solution() const;

            //- Return const access to the output properties dictionary
            inline const IOdictionary& outputProperties() const;

            //- Return output properties dictionary
            inline IOdictionary& outputProperties();

            //- Return the acceleration due to gravity
            inline const dimensionedVector& g() const;


            // Film geometry

                //- Return the patch normal vectors
                const volVectorField& nHat() const;

                //- Return the face area magnitudes [m^2]
                const volScalarField::Internal& magSf() const;

                //- Return the cell layer volume/area [m]
                const volScalarField& VbyA() const;

            // Addressing

                //- Return true if patchi on the film is a coupled
                //  patch to the primary film
                inline bool isCoupledPatch(const label filmPatchi) const;

                //- Return true if patchi on the primary film is a coupled
                //  patch to the film
                inline bool isFilmPatch(const label primaryPatchi) const;

                //- Return the list of patch IDs on the primary film coupled
                //  to this film
                inline const labelList& primaryPatchIDs() const;

                //- Return the list of patch IDs internally coupled with the
                //  primary film
                inline const labelList& intCoupledPatchIDs() const;

                //- Return film ID corresponding to primaryPatchID
                inline label filmPatchID(const label primaryPatchID) const;

                //- Return the list of patch IDs opposite to internally
                //  coupled patches
                const labelList& passivePatchIDs() const;


        // Patch type information

            //- Return boundary types for mapped field patches
            //  Also maps internal field value
            //  Mapping film prescribed by underlying mapped poly patch
            template<class Type>
            wordList mappedFieldAndInternalPatchTypes() const;

        // Helper

            //- Return the coupled patch ID paired with coupled patch
            //  filmPatchi
            label nbrCoupledPatchID
            (
                const surfaceFilm& nbrFilm,
                const label filmPatchi
            ) const;

            //- Convert a film field to the primary film
            template<class Type>
            void toPrimary
            (
                const label filmPatchi,
                Field<Type>& filmField
            ) const;

            //- Convert a primary film field to the film
            template<class Type>
            tmp<Field<Type>> toFilm
            (
                const label filmPatchi,
                const Field<Type>& primaryPatchField
            ) const;

            //- Convert a primary film field to the film
            template<class Type>
            tmp<Field<Type>> toFilm
            (
                const label filmPatchi,
                const tmp<Field<Type>> & primaryPatchField
            ) const;

            //- Return a primary patch field mapped the film
            template<class Type>
            void toFilm
            (
                Field<Type>& rf,
                const typename VolField<Type>::Boundary& pBf
            ) const;


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


        // Transfer fields - to the primary film

            //- Return mass transfer source - Eulerian phase only
            virtual tmp<volScalarField> primaryMassTrans() const = 0;

            //- Return the film mass available for transfer
            virtual const volScalarField& cloudMassTrans() const = 0;

            //- Return the parcel diameters originating from film
            virtual const volScalarField& cloudDiameterTrans() const = 0;


        // Evolution

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

            //- Main driver routing to evolve the film - calls other evolves
            virtual void evolve();

            //- Pre-evolve film
            virtual void preEvolveFilm();

            //- Evolve the film
            virtual void evolveFilm();

            //- Post-evolve film
            virtual void postEvolveFilm();


        // I-O

            //- Provide some feedback
            virtual void info() = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const surfaceFilm&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "surfaceFilmI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "surfaceFilmTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
