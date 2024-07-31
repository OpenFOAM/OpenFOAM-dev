/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2024 OpenFOAM Foundation
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
    Foam::fv::sixDoFAcceleration

Description
    Solid-body 6-DoF acceleration source

Usage
    Example usage:
    \verbatim
    SBM
    {
        type            sixDoFAcceleration;

        accelerations   <function1>;
    }
    \endverbatim

SourceFiles
    sixDoFAcceleration.C

\*---------------------------------------------------------------------------*/

#ifndef sixDoFAcceleration_H
#define sixDoFAcceleration_H

#include "fvModel.H"
#include "Function1.H"
#include "dimensionedTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                 Class sixDoFAcceleration Declaration
\*---------------------------------------------------------------------------*/

class sixDoFAcceleration
:
    public fvModel
{
private:

    // Private data

        //- Velocity field name, default = U
        word UName_;

        //- Linear acceleration function
        autoPtr<Function1<vector>> acceleration_;

        //- Angular velocity function
        autoPtr<Function1<vector>> angularVelocity_;

        //- Angular acceleration function
        autoPtr<Function1<vector>> angularAcceleration_;

        //- Optional gravitational acceleration
        dimensionedVector g_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs(const dictionary& dict);

        //- Add force to a momentum equation
        template<class AlphaFieldType, class RhoFieldType>
        inline void addForce
        (
            const AlphaFieldType& alpha,
            const RhoFieldType& rho,
            const volVectorField& U,
            fvMatrix<vector>& eqn
        ) const;


public:

    //- Runtime type information
    TypeName("sixDoFAcceleration");


    // Constructors

        //- Construct from components
        sixDoFAcceleration
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        sixDoFAcceleration
        (
            const sixDoFAcceleration&
        ) = delete;


    //- Destructor
    virtual ~sixDoFAcceleration()
    {}


    // Member Functions

        //- Return the list of fields for which the fvModel adds source term
        //  to the transport equation
        virtual wordList addSupFields() const;

        //- Source term to momentum equation
        virtual void addSup
        (
            const volVectorField& U,
            fvMatrix<vector>& eqn
        ) const;

        //- Source term to compressible momentum equation
        virtual void addSup
        (
            const volScalarField& rho,
            const volVectorField& U,
            fvMatrix<vector>& eqn
        ) const;

        //- Source term to phase momentum equation
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


        //- Read dictionary
        virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const sixDoFAcceleration&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
