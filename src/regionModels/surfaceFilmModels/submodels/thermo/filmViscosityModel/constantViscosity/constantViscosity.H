/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::constantViscosity

Description
    Constant viscosity model

SourceFiles
    constantViscosity.C

\*---------------------------------------------------------------------------*/

#ifndef constantViscosity_H
#define constantViscosity_H

#include "filmViscosityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                      Class constantViscosity Declaration
\*---------------------------------------------------------------------------*/

class constantViscosity
:
    public filmViscosityModel
{
private:

    // Private member functions

        //- Disallow default bitwise copy construct
        constantViscosity(const constantViscosity&);

        //- Disallow default bitwise assignment
        void operator=(const constantViscosity&);


protected:

    // Protected data

        //- Constant viscosity [Pa.s]
        scalar mu0_;



public:

    //- Runtime type information
    TypeName("constant");


    // Constructors

        //- Construct from surface film model
        constantViscosity
        (
            surfaceFilmModel& owner,
            const dictionary& dict,
            volScalarField& mu
        );


    //- Destructor
    virtual ~constantViscosity();


    // Member Functions

        // Evolution

            //- Correct
            virtual void correct
            (
                const volScalarField& p,
                const volScalarField& T
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
