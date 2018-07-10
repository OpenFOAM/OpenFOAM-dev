/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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
    Foam::ParticleStressModels::HarrisCrighton

Description
    Inter-particle stress model of Harris and Crighton

    The stress value takes the following form:
    \f[
        \frac{P_s \alpha^\beta}{ \mathrm{max} \left( \alpha_{pack} - \alpha ,
        \epsilon ( 1 - \alpha ) \right) }
    \f]
    Here, \f$\alpha\f$ is the volume fraction of the dispersed phase, and the
    other values are modelling constants. A small value \f$\epsilon\f$ is used
    to limit the denominator to ensure numerical stability.

    Reference:
    \verbatim
        "Solitons, solitary waves, and voidage disturbances in gas-fluidized
        beds"
        S Harris and D Crighton,
        Journal of Fluid Mechanics
        Volume 266, Pages 243-276, 1994
    \endverbatim

SourceFiles
    HarrisCrighton.C

\*---------------------------------------------------------------------------*/

#ifndef HarrisCrighton_H
#define HarrisCrighton_H

#include "ParticleStressModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace ParticleStressModels
{

/*---------------------------------------------------------------------------*\
                        Class HarrisCrighton Declaration
\*---------------------------------------------------------------------------*/

class HarrisCrighton
:
    public ParticleStressModel
{
    // Private data

        //- Solid pressure coefficient
        scalar pSolid_;

        //- Exponent of the volume fraction
        scalar beta_;

        //- Smallest allowable difference from the packed volume fraction
        scalar eps_;


    // Private member functions

        //- Return the limited denominator of the radial distribution function
        tmp<Field<scalar>> denominator(const Field<scalar>& alpha) const;


public:

    //- Runtime type information
    TypeName("HarrisCrighton");


    //- Constructors

        //- Construct from components
        HarrisCrighton(const dictionary& dict);

        //- Construct copy
        HarrisCrighton(const HarrisCrighton& hc);

        //- Clone
        virtual autoPtr<ParticleStressModel> clone() const
        {
            return autoPtr<ParticleStressModel>
            (
                new HarrisCrighton(*this)
            );
        }


    //- Destructor
    virtual ~HarrisCrighton();


    //- Member Functions

        //- Collision stress
        tmp<Field<scalar>> tau
        (
            const Field<scalar>& alpha,
            const Field<scalar>& rho,
            const Field<scalar>& uRms
        ) const;

        //- Collision stress derivaive w.r.t. the volume fraction
        tmp<Field<scalar>> dTaudTheta
        (
            const Field<scalar>& alpha,
            const Field<scalar>& rho,
            const Field<scalar>& uRms
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace ParticleStressModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
