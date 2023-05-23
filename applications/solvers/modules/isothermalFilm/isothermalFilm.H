/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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
    Foam::solvers::isothermalFilm

Description
    Solver module for flow of compressible isothermal liquid films

    Uses the flexible PIMPLE (PISO-SIMPLE) solution for time-resolved and
    pseudo-transient and steady simulations.

    Optional fvModels and fvConstraints are provided to enhance the simulation
    in many ways including adding various sources, Lagrangian
    particles, surface film etc. and constraining or limiting the solution.

SourceFiles
    isothermalFilm.C

See also
    Foam::solver

\*---------------------------------------------------------------------------*/

#ifndef isothermalFilm_H
#define isothermalFilm_H

#include "solver.H"
#include "rhoThermo.H"
#include "filmCompressibleMomentumTransportModel.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class surfaceTensionModel;
class mappedPatchBase;

namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class isothermalFilm Declaration
\*---------------------------------------------------------------------------*/

class isothermalFilm
:
    public solver
{
    // Control parameters

        //- Maximum allowed Courant number
        scalar maxCo;

        //- Maximum time-step
        scalar maxDeltaT_;


    // Continuity properties

        //- Current maximum Courant number for time-step control
        scalar CoNum;

        //- Current cumulative continuity error
        scalar cumulativeContErr;


protected:

    // Thermophysical properties

        //- Pointer to the fluid thermophysical properties
        autoPtr<rhoThermo> thermoPtr_;

        //- Reference to the fluid thermophysical properties
        rhoThermo& thermo_;

        //- The thermodynamic pressure field
        volScalarField& p;


    // Film

        //- List of film wall patch IDs
        labelList wallPatchIDs;

        //- Film surface patch ID
        //  Set to -1 if the surface patch is empty or not coupled
        label surfacePatchID;

        //- Film wall normal
        volVectorField nHat_;

        //- Film cell cross-sectional area magnitude
        volScalarField::Internal magSf_;

        //- Film cell volume/wall face area
        volScalarField VbyA_;

        //- Bool returned by initFilmMesh()
        bool initialised_;

        //- Film thickness
        volScalarField delta_;

        //- Film volume fraction in the cell layer
        volScalarField alpha_;

        //- Film thickness below which the surface is considered dry
        dimensionedScalar deltaWet;


    // Kinematic properties

        //- Film velocity field
        volVectorField U_;

        //- Film mass-flux field
        surfaceScalarField alphaRhoPhi_;

        //- Film volumetric-flux field
        surfaceScalarField phi_;


    // Interface properties

        //- Pointer to the surface tension coefficient model
        autoPtr<surfaceTensionModel> surfaceTension;

        //- Set true if the surface tension coefficient is non-uniform
        //  to include the thermocapillary force
        bool thermocapillary;


    // Momentum transport

        //- Pointer to the momentum transport model
        autoPtr<filmCompressible::momentumTransportModel> momentumTransport;


    // Cached temporary fields

        //- Cached momentum matrix
        //  shared between the momentum predictor and pressure corrector
        tmp<fvVectorMatrix> tUEqn;

        //- Continuity error
        tmp<volScalarField::Internal> contErr;


private:

    // Private Member Functions

        //- Initialise film specific mesh topology and geometry
        bool initFilmMesh();

        //- Helper function to map the delta BC types
        //  to the corresponding alpha BC types
        wordList alphaTypes() const;

        //- Constrain the given field such that the flow
        //  is in the plane of the film.
        //  This involves setting the flux on the top and bottom patches to 0
        template<class FieldType>
        void constrainField(FieldType& field) const;

        //- Constrain a field such that the flow
        //  is in the plane of the film.
        //  This involves setting the flux on the top and bottom patches to 0
        template<class FieldType>
        tmp<FieldType> constrainedField(const FieldType& field) const;

        //- Constrain a field such that the flow
        //  is in the plane of the film.
        //  This involves setting the flux on the top and bottom patches to 0
        template<class FieldType>
        tmp<FieldType> constrainedField(const tmp<FieldType>& tfield) const;

        //- Correct the cached Courant number
        void correctCoNum();

        //  Solve the explicit continuity equation for the film volume-fraction
        //  to predict the film thickness
        void continuityPredictor();

        // Calculate the continuity error caused by limiting alpha
        void correctContinuityError();

        //- Print the continuity errors
        void continuityErrors();

        //- Update film thickness delta from the film volume-fraction
        void correctDelta();

        //- Construct the film volume-fraction elliptic equation
        //  and correct the film thickness
        void correctAlpha();

        //- Buoyant pressure divided by alpha*density
        tmp<surfaceScalarField> pbByAlphaRhof() const;

        //- Buoyant pressure divided by alpha
        tmp<surfaceScalarField> pbByAlphaf() const;

        //- Buoyant pressure divided by alpha*grad(rho)
        tmp<surfaceScalarField> pbByAlphaGradRhof() const;

        //- Capillary pressure
        tmp<volScalarField> pc(const volScalarField& sigma) const;

        //- External pressure
        //    Adjacent region fluid pressure
        //    Accumulated Lagrangian particle impingement pressure
        tmp<volScalarField> pe() const;

        //- Contact force obtained from the given surface tension coefficient
        //  and the contact angle distribution from the wall delta patch field
        tmp<volVectorField::Internal> contactForce
        (
            const volScalarField& sigma
        ) const;


protected:

    // Protected Member Functions

        //- Read controls
        void readControls();


public:

    // Public Data

        //- Acceleration due to gravity
        const uniformDimensionedVectorField g;

        //- Film wall normal
        const volVectorField& nHat;

        //- Film cell cross-sectional area magnitude
        const volScalarField::Internal& magSf;

        //- Film cell volume/wall face area
        const volScalarField& VbyA;

        //- Film thickness
        const volScalarField& delta;

        //- Film volume fraction in the cell layer
        const volScalarField& alpha;

        //- Reference to the fluid thermophysical properties
        const rhoThermo& thermo;

        //- Reference to the thermodynamic density field
        const volScalarField& rho;

        //- Reference to the film velocity field
        const volVectorField& U;

        //- Reference to the film mass-flux field
        const surfaceScalarField& alphaRhoPhi;

        //- Reference to the film volumetric-flux field
        const surfaceScalarField& phi;


public:

    //- Runtime type information
    TypeName("isothermalFilm");


    // Constructors

        //- Construct from region mesh and thermophysical properties
        isothermalFilm(fvMesh& mesh, autoPtr<rhoThermo>);

        //- Construct from region mesh
        isothermalFilm(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        isothermalFilm(const isothermalFilm&) = delete;


    //- Destructor
    virtual ~isothermalFilm();


    // Member Functions

        //- Return the film surface patch
        const fvPatch& surfacePatch() const;

        //- Return the film surface patch region-region map
        const mappedPatchBase& surfacePatchMap() const;

        //- Return the film surface tension coefficient field
        tmp<volScalarField> sigma() const;

        //- Return the current maximum time-step for stable solution
        virtual scalar maxDeltaT() const;

        //- Called at the start of the time-step, before the PIMPLE loop
        virtual void preSolve();

        //- Called at the start of the PIMPLE loop to move the mesh
        virtual void moveMesh();

        //- Called at the start of the PIMPLE loop
        virtual void prePredictor();

        //- Construct and optionally solve the momentum equation
        virtual void momentumPredictor();

        //- Construct and solve the energy equation,
        //  convert to temperature
        //  and update thermophysical and transport properties
        virtual void thermophysicalPredictor();

        //- Construct and solve the pressure equation in the PISO loop
        virtual void pressureCorrector();

        //- Correct the momentum and thermophysical transport modelling
        virtual void postCorrector();

        //- Called after the PIMPLE loop at the end of the time-step
        virtual void postSolve();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const isothermalFilm&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "isothermalFilmTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
