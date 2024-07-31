/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2024 OpenFOAM Foundation
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
    Foam::fv::OUForce

Description
    Calculates and applies the random OU (Ornstein-Uhlenbeck) process force to
    the momentum equation for direct numerical simulation of boxes of isotropic
    turbulence.

    The energy spectrum is calculated and written at write-times which is
    particularly useful to test and compare LES SGS models.

Note
    This random OU process force uses a FFT to generate the force field which
    is not currently parallelised.  Also the mesh the FFT is applied to must
    be isotropic and have a power of 2 cells in each direction.

Usage
    Example usage:
    \verbatim
    OUForce
    {
        type    OUForce;

        libs    ("librandomProcesses.so");

        sigma   0.090295;
        alpha   0.81532;
        kUpper  10;
        kLower  7;
    }
    \endverbatim

SourceFiles
    OUForce.C

\*---------------------------------------------------------------------------*/

#ifndef OUForce_H
#define OUForce_H

#include "fvModel.H"
#include "Kmesh.H"
#include "OUprocess.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class OUForce Declaration
\*---------------------------------------------------------------------------*/

class OUForce
:
    public fvModel
{
    // Private Data

        //- Name of the velocity field
        word UName_;

        //- The wavenumber mesh used for the random OU process
        Kmesh K_;

        //- The random OU process force generator
        OUprocess forceGen_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs(const dictionary& dict);


public:

    //- Runtime type information
    TypeName("OUForce");


    // Constructors

        //- Construct from explicit source name and mesh
        OUForce
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        OUForce(const OUForce&) = delete;


    // Member Functions

        // Checks

            //- Return the list of fields for which the fvModel adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;


        // Evaluate

            //- Add explicit contribution to incompressible momentum equation
            virtual void addSup
            (
                const volVectorField& field,
                fvMatrix<vector>& eqn
            ) const;


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

            //- Write energy spectrum
            virtual bool write(const bool write = true) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const OUForce&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
