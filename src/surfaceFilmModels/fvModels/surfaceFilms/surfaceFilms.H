/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2022 OpenFOAM Foundation
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
    Foam::fv::surfaceFilms

Description
    Surface films fvModel

Usage
    Example usage:
    \verbatim
    surfaceFilms
    {
        type    surfaceFilms;
        libs    ("libsurfaceFilmModels.so");
    }
    \endverbatim

SourceFiles
    surfaceFilms.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceFilms_H
#define surfaceFilms_H

#include "fvModel.H"
#include "thermoSurfaceFilm.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                          Class surfaceFilms Declaration
\*---------------------------------------------------------------------------*/

class surfaceFilms
:
    public fvModel
{
    // Private Data

        //- The surface film models
        PtrList<regionModels::thermoSurfaceFilm> surfaceFilms_;

        //- For each surface film model, the corresponding density name in the
        //  primary region
        wordList surfaceFilmPrimaryRhoNames_;

        //- List of fields for which the fvModel adds source term
        //  to the transport equation
        wordList fieldNames_;

        //- Current time index (used for updating)
        mutable label curTimeIndex_;


public:

    //- Runtime type information
    TypeName("surfaceFilms");


    // Constructors

        //- Construct from explicit source name and mesh
        surfaceFilms
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construction
        surfaceFilms
        (
            const surfaceFilms&
        ) = delete;


    // Member Functions

        // Checks

            //- Return the list of fields for which the option adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;

            //- Return the maximum time-step for stable operation
            virtual scalar maxDeltaT() const;


        // Correct

            //- Solve the Lagrangian surfaceFilms and update the sources
            virtual void correct();


        // Add explicit and implicit contributions to compressible equation

            //- Add source to continuity equation
            virtual void addSup
            (
                fvMatrix<scalar>& eqn,
                const word& fieldName
            ) const;

            //- Add source to pressure or enthalpy equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const word& fieldName
            ) const;

            //- Add source to momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const word& fieldName
            ) const;


        // Mesh changes

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&);

            //- Update for mesh motion
            virtual bool movePoints();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const surfaceFilms&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
