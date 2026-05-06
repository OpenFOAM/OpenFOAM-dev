/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::fv::solidElectricalConduction

Description
    This fvModel represents conduction through a solid region. It solves for
    the electric potential, calculates the current, and applies the
    corresponding Joule heating source to the energy equation. It requires the
    user to specify the electric potential field \c phi in the zero directory.

Usage
    Example usage with a constant scalar conductivity:
    \verbatim
    solidElectricalConduction
    {
        type            solidElectricalConduction;

        sigma           35.5e6; // <-- 'sigma<scalar>' is also accepted
    }
    \endverbatim

    Example usage with a constant tensor conductivity:
    \verbatim
    solidElectricalConduction
    {
        type            solidElectricalConduction;

        sigma<tensor>   (35.5e6 0 0 0 35.5e6 0 0 0 35.5e6);
    }
    \endverbatim

    Example usage with a functional scalar conductivity, and writing the
    conductivity out for post-processing:
    \verbatim
    solidElectricalConduction
    {
        type            solidElectricalConduction;

        sigma
        {
            type        distanceFunction;

            function
            {
                type        polynomial;
                coeffs      (35.5e6 1e9);
            }

            direction   (0 1 0);
        }

        writeSigma      yes;
    }
    \endverbatim

SourceFiles
    solidElectricalConduction.C

\*---------------------------------------------------------------------------*/

#ifndef solidElectricalConduction_H
#define solidElectricalConduction_H

#include "fvModel.H"
#include "volFields.H"
#include "FunctionalGeometricFieldFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class solidElectricalConduction Declaration
\*---------------------------------------------------------------------------*/

class solidElectricalConduction
:
    public fvModel
{
    // Private Data

        //- The electric potential
        volScalarField phi_;

        //- The current
        surfaceScalarField I_;

        //- The scalar conductivity
        autoPtr<FunctionalGeometricField<scalar, fvMesh>> sigmaScalarPtr_;

        //- The tensor conductivity
        autoPtr<FunctionalGeometricField<tensor, fvMesh>> sigmaTensorPtr_;

        //- Write the conductivity?
        bool writeSigma_;

        //- Has the mesh changed since the last solution?
        bool meshChanged_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs(const dictionary& dict);


public:

    //- Runtime type information
    TypeName("solidElectricalConduction");


    // Constructors

        //- Construct from explicit source name and mesh
        solidElectricalConduction
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        solidElectricalConduction(const solidElectricalConduction&) = delete;


    // Member Functions

        // Checks

            //- Return the list of fields for which the fvModel adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;


        // Evaluate

            //- Add explicit contribution to compressible energy equation
            virtual void addSup
            (
                const volScalarField& rho,
                const volScalarField& he,
                fvMatrix<scalar>& eqn
            ) const;


        //- Solve for the electric potential
        virtual void correct();


        // Mesh changes

            //- Update for mesh motion
            virtual bool movePoints();

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&);


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);

            //- Write fvModel data
            virtual bool write(const bool write = true) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const solidElectricalConduction&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
