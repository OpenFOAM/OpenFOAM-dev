/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2021 OpenFOAM Foundation
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
    Foam::fixedTrim

Description
    Fixed trim coefficients

SourceFiles
    fixedTrim.C

\*---------------------------------------------------------------------------*/

#ifndef fixedTrim_H
#define fixedTrim_H

#include "trimModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class fixedTrim Declaration
\*---------------------------------------------------------------------------*/

class fixedTrim
:
    public trimModel
{

protected:

    // Protected data

        //- Geometric angle of attack [rad]
        scalarField thetag_;


public:

    //- Run-time type information
    TypeName("fixedTrim");

    //- Constructor
    fixedTrim(const fv::rotorDiskSource& rotor, const dictionary& dict);

    //- Destructor
    virtual ~fixedTrim();


    // Member Functions

        //- Read
        void read(const dictionary& dict);

        //- Return the geometric angle of attack [rad]
        virtual tmp<scalarField> thetag() const;

        //- Correct the model
        virtual void correct
        (
            const vectorField& U,
            vectorField& force
        ) const;

        //- Correct the model for compressible flow
        virtual void correct
        (
            const volScalarField rho,
            const vectorField& U,
            vectorField& force
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
