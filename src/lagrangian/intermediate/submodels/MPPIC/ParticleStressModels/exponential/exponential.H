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
    Foam::ParticleStressModels::exponential

Description
    Exponential inter-particle stress model of the same form as used in
    twoPhaseEulerFoam

SourceFiles
    exponential.C

\*---------------------------------------------------------------------------*/

#ifndef exponential_H
#define exponential_H

#include "ParticleStressModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace ParticleStressModels
{

/*---------------------------------------------------------------------------*\
                        Class exponential Declaration
\*---------------------------------------------------------------------------*/

class exponential
:
    public ParticleStressModel
{
    // Private data

        //- Pre-exponential factor
        scalar preExp_;

        //- Maximum limit of the exponential
        scalar expMax_;

        //- Front coefficient
        scalar g0_;


public:

    //- Runtime type information
    TypeName("exponential");


    //- Constructors

        //- Construct from components
        exponential(const dictionary& dict);

        //- Construct copy
        exponential(const exponential& hc);

        //- Clone
        virtual autoPtr<ParticleStressModel> clone() const
        {
            return autoPtr<ParticleStressModel>
            (
                new exponential(*this)
            );
        }


    //- Destructor
    virtual ~exponential();


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
