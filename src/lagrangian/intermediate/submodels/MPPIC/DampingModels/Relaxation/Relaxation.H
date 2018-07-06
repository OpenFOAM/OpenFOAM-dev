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
    Foam::DampingModels::Relaxation

Description
    Relaxation collisional damping model.

    Particle velocities are relaxed towards the local mean over a time-scale.

    Reference:
    \verbatim
        "An improved collision damping time for MP-PIC calculations of dense
        particle flows with applications to polydisperse sedimenting beds and
        colliding particle jets"
        P O'Rourke and D Snider
        Chemical Engineering Science
        Volume 65, Issue 22, Pages 6014-6028, November 2010
    \endverbatim

SourceFiles
    Relaxation.C

\*---------------------------------------------------------------------------*/

#ifndef Relaxation_H
#define Relaxation_H

#include "DampingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace DampingModels
{

/*---------------------------------------------------------------------------*\
                         Class Relaxation Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class Relaxation
:
    public DampingModel<CloudType>
{
private:

    // Private data

        //- Velocity average
        const AveragingMethod<vector>* uAverage_;

        //- Reciprocal of the time scale average
        autoPtr<AveragingMethod<scalar>> oneByTimeScaleAverage_;


public:

    //- Runtime type information
    TypeName("relaxation");

    // Constructors

        //- Construct from components
        Relaxation(const dictionary& dict, CloudType& owner);

        //- Construct copy
        Relaxation(const Relaxation<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<DampingModel<CloudType>> clone() const
        {
            return autoPtr<DampingModel<CloudType>>
            (
                new Relaxation<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~Relaxation();


    //- Member Functions

        //- Calculate the damping time scales
        virtual void cacheFields(const bool store);

        //- Calculate the velocity correction
        virtual vector velocityCorrection
        (
            typename CloudType::parcelType& p,
            const scalar deltaT
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace DampingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Relaxation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
