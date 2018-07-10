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
    Foam::IsotropyModels::Stochastic

Description
    Stochastic return-to-isotropy model.

    Particle velocities are modified by sampling a gaussian-plus-delta
    distribution, which depends on a time-scale. This randomises some particle
    velocities whilst leaving others unchanged. The lower the value of the
    time-scale, the greater the proportion of particle velocities affected.

    A correction step is performed at the end to ensure that the model
    conserves both momentum and granular temperature.

    Reference:
    \verbatim
        "Inclusion of collisional return-to-isotropy in the MP-PIC method"
        P O'Rourke and D Snider
        Chemical Engineering Science
        Volume 80, Issue 0, Pages 39-54, December 2012
    \endverbatim

SourceFiles
    Stochastic.C

\*---------------------------------------------------------------------------*/

#ifndef Stochastic_H
#define Stochastic_H

#include "IsotropyModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace IsotropyModels
{

/*---------------------------------------------------------------------------*\
                         Class Stochastic Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class Stochastic
:
    public IsotropyModel<CloudType>
{
private:

    // Private Member Functions

        //- Sample a gaussian distribution using the Box-Muller method
        scalar sampleGauss();


public:

    //- Runtime type information
    TypeName("stochastic");

    // Constructors

        //- Construct from components
        Stochastic(const dictionary& dict, CloudType& owner);

        //- Construct copy
        Stochastic(const Stochastic<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<IsotropyModel<CloudType>> clone() const
        {
            return autoPtr<IsotropyModel<CloudType>>
            (
                new Stochastic<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~Stochastic();


    //- Member Functions

        //- Calculate velocities
        virtual void calculate();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace IsotropyModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Stochastic.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
