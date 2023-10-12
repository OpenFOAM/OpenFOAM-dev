/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2023 OpenFOAM Foundation
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
    Foam::fv::interRegionPorosityForce

Description
    This model applies the force exerted on the fluid by a porous media, the
    extent of which is defined by an overlapping region

Usage
    Example usage, here employing the Darcy-Forchheimer model:
    \verbatim
    interRegionPorosityForceCoeffs
    {
        type            DarcyForchheimer;

        DarcyForchheimerCoeffs
        {
            d   d [0 -2 0 0 0 0 0] (5e7 -1000 -1000);
            f   f [0 -1 0 0 0 0 0] (0 0 0);

            coordinateSystem
            {
                type    cartesian;
                origin  (0 0 0);
                coordinateRotation
                {
                    type    axesRotation;
                    e1  (0.70710678 0.70710678 0);
                    e2  (0 0 1);
                }
            }
        }
    }
    \endverbatim

SourceFiles
    interRegionPorosityForce.C

\*---------------------------------------------------------------------------*/

#ifndef interRegionPorosityForce_H
#define interRegionPorosityForce_H

#include "interRegionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class porosityModel;

namespace fv
{


/*---------------------------------------------------------------------------*\
              Class interRegionPorosityForce Declaration
\*---------------------------------------------------------------------------*/

class interRegionPorosityForce
:
    public interRegionModel
{
    // Private data

        //- Velocity field name, default = U
        word UName_;

        //- Porosity mapping filter
        volScalarField::Internal filter_;

        //- Run-time selectable porosity model
        mutable autoPtr<porosityModel> porosityPtr_;


    // Protected Member Functions

        //- Non-virtual read
        void readCoeffs();


public:

    //- Runtime type information
    TypeName("interRegionPorosityForce");


    // Constructors

        //- Construct from components
        interRegionPorosityForce
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        interRegionPorosityForce
        (
            const interRegionPorosityForce&
        ) = delete;


    //- Destructor
    virtual ~interRegionPorosityForce()
    {}


    // Member Functions

        // Checks

            //- Return the list of fields for which the fvModel adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;


        // Evaluation

            //- Add source to momentum equation
            virtual void addSup
            (
                const volVectorField& U,
                fvMatrix<vector>& eqn
            ) const;

            //- Add source to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                const volVectorField& U,
                fvMatrix<vector>& eqn
            ) const;

            //- Add source to phase momentum equation
            virtual void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                const volVectorField& U,
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

            //- Read dictionary
            virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const interRegionPorosityForce&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
