/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
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
    Foam::waveModels::Airy

Description
    First-order wave model.

    Reference:
    \verbatim
        "On the Theory of Oscillatory Waves"
        G G Stokes
        Transactions of the Cambridge Philosophical Society (1847), Volume 8,
        Pages 441-455, Equations 18 and 19 (leading terms only)
    \endverbatim

SourceFiles
    Airy.C

\*---------------------------------------------------------------------------*/

#ifndef Airy_H
#define Airy_H

#include "waveModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveModels
{

/*---------------------------------------------------------------------------*\
                           Class Airy Declaration
\*---------------------------------------------------------------------------*/

class Airy
:
    public waveModel
{
public:

    //- Runtime type information
    TypeName("Airy");


    // Constructors

        //- Construct from a database and a dictionary
        Airy(const objectRegistry& db, const dictionary& dict);

        //- Construct a clone
        virtual autoPtr<waveModel> clone() const
        {
            return autoPtr<waveModel>(new Airy(*this));
        }


    //- Destructor
    virtual ~Airy();


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

} // End namespace waveModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
