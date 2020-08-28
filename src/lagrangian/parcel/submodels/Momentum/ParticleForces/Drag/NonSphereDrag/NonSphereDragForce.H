/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::NonSphereDragForce

Description
    Drag model for non-spherical particles.

    The drag coefficient is given by:
    \f[
        Cd = \frac{24}{Re} (1 + a Re^{b}) + Re \frac{c}{Re + d}
    \f]

    Where \f$a\f$, \f$b\f$, \f$c\f$, and \f$d\f$ are coefficients, calculated
    as functions of \f$phi\f$. \f$phi\f$ is the ratio of the surface area of a
    sphere with the same volume as the particle to the actual surface area of
    the particle, and must be between 0 and 1.

    Reference:
    \verbatim
        Haider, A., & Levenspiel, O. (1989).
        Drag coefficient and terminal velocity of spherical and nonspherical
        particles.
        Powder technology, 58(1), 63-70.
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef NonSphereDragForce_H
#define NonSphereDragForce_H

#include "ParticleForce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class NonSphereDragForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NonSphereDragForce
:
    public ParticleForce<CloudType>
{
protected:

    // Protected Data

        //- The ratio of the surface area of a sphere with the same volume as
        //  the particle to the actual surface area of the particle. Between 0
        //  and 1.
        const scalar phi_;

        //- Coefficient
        const scalar a_;

        //- Coefficient
        const scalar b_;

        //- Coefficient
        const scalar c_;

        //- Coefficient
        const scalar d_;


public:

    //- Runtime type information
    TypeName("nonSphereDrag");


    // Constructors

        //- Construct from mesh
        NonSphereDragForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Construct copy
        NonSphereDragForce(const NonSphereDragForce<CloudType>& df);

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType>> clone() const
        {
            return autoPtr<ParticleForce<CloudType>>
            (
                new NonSphereDragForce<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~NonSphereDragForce();


    // Member Functions

        // Evaluation

            //- Calculate the non-coupled force
            virtual forceSuSp calcCoupled
            (
                const typename CloudType::parcelType& p,
                const typename CloudType::parcelType::trackingData& td,
                const scalar dt,
                const scalar mass,
                const scalar Re,
                const scalar muc
            ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "NonSphereDragForce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
