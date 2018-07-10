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
    Foam::PackingModels::Explicit

Description
    Explicit model for applying an inter-particle stress to the particles.

    The inter-particle stress is calculated using current particle locations.
    This force is then applied only to the particles that are moving towards
    regions of close pack. The resulting velocity change is limited using an
    abtracted correction velocity limiter.

    Reference:
    \verbatim
        "An Incompressible Three-Dimensional Multiphase Particle-in-Cell Model
        for Dense Particle Flows"
        D Snider
        Journal of Computational Physics
        Volume 170, Issue 2, Pages 523-549, July 2001
    \endverbatim

SourceFiles
    Explicit.C

\*---------------------------------------------------------------------------*/

#ifndef Explicit_H
#define Explicit_H

#include "PackingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PackingModels
{

/*---------------------------------------------------------------------------*\
                         Class Explicit Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class Explicit
:
    public PackingModel<CloudType>
{
private:

    //- Private data

        //- Volume fraction average
        const AveragingMethod<scalar>* volumeAverage_;

        //- Velocity average
        const AveragingMethod<vector>* uAverage_;

        //- Stress average field
        autoPtr<AveragingMethod<scalar>> stressAverage_;

        //- Correction limiter
        autoPtr<CorrectionLimitingMethod> correctionLimiting_;


public:

    //- Runtime type information
    TypeName("explicit");

    // Constructors

        //- Construct from components
        Explicit(const dictionary& dict, CloudType& owner);

        //- Construct copy
        Explicit(const Explicit<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<PackingModel<CloudType>> clone() const
        {
            return autoPtr<PackingModel<CloudType>>
            (
                new Explicit<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~Explicit();


    // Member Functions

        //- Calculate the inter particles stresses
        virtual void cacheFields(const bool store);

        //- Calculate the velocity correction
        virtual vector velocityCorrection
        (
            typename CloudType::parcelType& p,
            const scalar deltaT
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PackingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Explicit.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
