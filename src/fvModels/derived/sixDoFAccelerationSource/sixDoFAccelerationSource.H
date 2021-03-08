/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2021 OpenFOAM Foundation
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
    Foam::fv::sixDoFAccelerationSource

Description
    Solid-body 6-DoF acceleration source

Usage
    Example usage:
    \verbatim
    SBM
    {
        type            sixDoFAccelerationSource;

        accelerations   <function1>;
    }
    \endverbatim

SourceFiles
    sixDoFAccelerationSource.C

\*---------------------------------------------------------------------------*/

#ifndef sixDoFAccelerationSource_H
#define sixDoFAccelerationSource_H

#include "fvModel.H"
#include "Function1.H"
#include "dimensionedTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                 Class sixDoFAccelerationSource Declaration
\*---------------------------------------------------------------------------*/

class sixDoFAccelerationSource
:
    public fvModel
{
public:

        //- Accelerations type comprising
        //      linear acceleration vector
        //      angular velocity vector
        //      angular acceleration vector
        typedef Vector<vector> accelerationVectors;


private:

    // Private data

        //- Velocity field name, default = U
        word UName_;

        //- Accelerations function returning a vector containing
        //      linear acceleration vector
        //      angular velocity vector
        //      angular acceleration vector
        autoPtr<Function1<accelerationVectors>> accelerations_;

        //- Optional gravitational acceleration
        dimensionedVector g_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs();

        //- Source term to momentum equation
        template<class RhoFieldType>
        void addSup
        (
            const RhoFieldType& rho,
            fvMatrix<vector>& eqn,
            const word& fieldName
        ) const;


public:

    //- Runtime type information
    TypeName("sixDoFAccelerationSource");


    // Constructors

        //- Construct from components
        sixDoFAccelerationSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construction
        sixDoFAccelerationSource
        (
            const sixDoFAccelerationSource&
        ) = delete;


    //- Destructor
    virtual ~sixDoFAccelerationSource()
    {}


    // Member Functions

        //- Return the list of fields for which the fvModel adds source term
        //  to the transport equation
        virtual wordList addSupFields() const;

        //- Source term to momentum equation
        virtual void addSup
        (
            fvMatrix<vector>& eqn,
            const word& fieldName
        ) const;

        //- Source term to compressible momentum equation
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<vector>& eqn,
            const word& fieldName
        ) const;

        //- Read dictionary
        virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const sixDoFAccelerationSource&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "sixDoFAccelerationSourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
