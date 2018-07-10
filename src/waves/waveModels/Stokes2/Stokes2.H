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
    Foam::waveModels::Stokes2

Description
    Second-order wave model.

    Reference:
    \verbatim
        "On the Theory of Oscillatory Waves"
        G G Stokes
        Transactions of the Cambridge Philosophical Society (1847), Volume 8,
        Pages 441-455, Equations 18 and 19
    \endverbatim

SourceFiles
    Stokes2.C

\*---------------------------------------------------------------------------*/

#ifndef Stokes2_H
#define Stokes2_H

#include "Airy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveModels
{

/*---------------------------------------------------------------------------*\
                           Class Stokes2 Declaration
\*---------------------------------------------------------------------------*/

class Stokes2
:
    public Airy
{
public:

   //- Runtime type information
    TypeName("Stokes2");


    // Constructors

        //- Construct from a database and a dictionary
        Stokes2(const objectRegistry& db, const dictionary& dict);

        //- Construct a clone
        virtual autoPtr<waveModel> clone() const
        {
            return autoPtr<waveModel>(new Stokes2(*this));
        }


    //- Destructor
    virtual ~Stokes2();


    // Member Functions

        //- Get the wave elevation at a given time, mean velocity and local
        //  coordinates. Local x is aligned with the mean velocity.
        virtual tmp<scalarField> elevation
        (
            const scalar t,
            const scalar u,
            const scalarField& x
        ) const;

        //- Get the wave velocity at a given time, mean velocity and local
        //  coordinates. Local x is aligned with the mean velocity, and z with
        //  negative gravity.
        virtual tmp<vector2DField> velocity
        (
            const scalar t,
            const scalar u,
            const vector2DField& xz
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
} // End namespace waveModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
