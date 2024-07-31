/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2024 OpenFOAM Foundation
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
    Foam::fv::forcing

Description
    Base fvModel for forcing functions.

See also
    Foam::fv::isotropicDamping
    Foam::fv::verticalDamping

SourceFiles
    forcing.C

\*---------------------------------------------------------------------------*/

#ifndef forcing_H
#define forcing_H

#include "fvModel.H"
#include "Function1.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                           Class forcing Declaration
\*---------------------------------------------------------------------------*/

class forcing
:
    public fvModel
{
protected:

    // Protected Data

        //- Optionally write the forcing fields
        bool writeForceFields_;

        //- Forcing coefficient [1/s]
        mutable dimensionedScalar lambda_;

        //- Optional boundary forcing coefficient [1/s]
        mutable dimensionedScalar lambdaBoundary_;

        //- The scaling function
        autoPtr<Function1<scalar>> scale_;

        //- Origins of the scaling coordinate
        vectorField origins_;

        //- Directions of increasing scaling coordinate
        vectorField directions_;


    // Protected Member Functions

        //- Non-virtual read
        void readCoeffs(const dictionary& dict);

        //- Read the forcing coefficients
        void readLambda(const dictionary& dict);

        //- Calculate and return the volume average forcing region length
        dimensionedScalar regionLength() const;

        //- Return the scale distribution
        virtual tmp<volScalarField::Internal> scale() const;

        //- Return the force coefficient
        virtual tmp<volScalarField::Internal> forceCoeff() const;

        //- Optionally write the forcing fields:
        //  forcing:scale, forcing:forceCoeff
        void writeForceFields() const;


public:

    //- Runtime type information
    TypeName("forcing");


    // Constructors

        //- Construct from components
        forcing
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~forcing()
    {}


    // Member Functions

        //- Read dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
