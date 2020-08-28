/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
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
    Foam::PlessisMasliyahDragForce

Description
    Plessis-Masliyah drag model for spheres.

    Reference:
    \verbatim
        Du Plessis, J. P., & Masliyah, J. H. (1988).
        Mathematical modelling of flow through consolidated isotropic porous
        media.
        Transport in Porous Media, 3(2), 145-161.
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef PlessisMasliyahDragForce_H
#define PlessisMasliyahDragForce_H

#include "DenseDragForce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class PlessisMasliyahDragForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class PlessisMasliyahDragForce
:
    public DenseDragForce<CloudType>
{
public:

    //- Runtime type information
    TypeName("PlessisMasliyahDrag");


    // Constructors

        //- Construct from mesh
        PlessisMasliyahDragForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Construct copy
        PlessisMasliyahDragForce(const PlessisMasliyahDragForce<CloudType>& df);

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType>> clone() const
        {
            return autoPtr<ParticleForce<CloudType>>
            (
                new PlessisMasliyahDragForce<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~PlessisMasliyahDragForce();


    // Member Functions

        // Evaluation

            //- Calculate the coupled force
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
    #include "PlessisMasliyahDragForce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
