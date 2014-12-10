/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2014 OpenFOAM Foundation
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
    Foam::ParticleStressModels::Lun

Description
    Inter-particle stress model of Lun et al

    The stress value takes the following form:
    \f[
        \left( \alpha \rho + \alpha^2 \rho (1 + e) \frac{3}{5}
        \left( 1 - \left( \frac{\alpha}{\alpha_{pack}} \right)^\frac{1}{3}
        \right) \right) \frac{1}{3} \sigma^2
    \f]
    Here, \f$\alpha\f$ is the volume fraction of the dispersed phase,
    \f$\rho\f$ is the density of the dispersed phase, \f$e\f$ is a coefficient
    of restitution, and \f$\sigma\f$ is the RMS velocity fluctuation.

    Reference:
    \verbatim
        "Kinetic theories for granular flow: inelastic particles in Couette
        flow and slightly inelastic particles in a general flowfield"
        C Lun, S Savage, G Jeffrey, N Chepurniy
        Journal of Fluid Mechanics
        Volume 140, Pages 223-256, 1984
    \endverbatim

SourceFiles
    Lun.C

\*---------------------------------------------------------------------------*/

#ifndef Lun_H
#define Lun_H

#include "ParticleStressModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace ParticleStressModels
{

/*---------------------------------------------------------------------------*\
                        Class Lun Declaration
\*---------------------------------------------------------------------------*/

class Lun
:
    public ParticleStressModel
{
    // Private data

        //- Coefficient of restitution
        scalar e_;

        //- Smallest allowable difference from the packed volume fraction
        scalar eps_;


public:

    //- Runtime type information
    TypeName("Lun");


    //- Constructors

        //- Construct from components
        Lun(const dictionary& dict);

        //- Construct copy
        Lun(const Lun& hc);

        //- Clone
        virtual autoPtr<ParticleStressModel> clone() const
        {
            return autoPtr<ParticleStressModel>
            (
                new Lun(*this)
            );
        }


    //- Destructor
    virtual ~Lun();


    //- Member Functions

        //- Collision stress
        tmp<Field<scalar> > tau
        (
            const Field<scalar>& alpha,
            const Field<scalar>& rho,
            const Field<scalar>& uRms
        ) const;

        //- Collision stress derivaive w.r.t. the volume fraction
        tmp<Field<scalar> > dTaudTheta
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
