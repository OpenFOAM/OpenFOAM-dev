/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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
    Foam::waveModels::Stokes5

Description
    Fifth-order wave model.

    Reference:
    \verbatim
        Fenton, J. D. (1985).
        A fifth-order Stokes theory for steady waves.
        Journal of waterway, port, coastal, and ocean engineering, 111(2),
        216-234.
    \endverbatim

SourceFiles
    Stokes5.C

\*---------------------------------------------------------------------------*/

#ifndef Stokes5_H
#define Stokes5_H

#include "Stokes2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveModels
{

/*---------------------------------------------------------------------------*\
                           Class Stokes5 Declaration
\*---------------------------------------------------------------------------*/

class Stokes5
:
    public Stokes2
{
public:

   //- Runtime type information
    TypeName("Stokes5");


    // Constructors

        //- Construct from a database and a dictionary
        Stokes5(const objectRegistry& db, const dictionary& dict);

        //- Construct a clone
        virtual autoPtr<waveModel> clone() const
        {
            return autoPtr<waveModel>(new Stokes5(*this));
        }


    //- Destructor
    virtual ~Stokes5();


    // Member Functions

        //- Get the wave elevation at a given time and local coordinates. Local
        //  x is aligned with the direction of propagation.
        virtual tmp<scalarField> elevation
        (
            const scalar t,
            const scalarField& x
        ) const;

        //- Get the wave velocity at a given time and local coordinates. Local
        //  x is aligned with the direction of propagation, and z with negative
        //  gravity.
        virtual tmp<vector2DField> velocity
        (
            const scalar t,
            const vector2DField& xz
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
} // End namespace waveModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
