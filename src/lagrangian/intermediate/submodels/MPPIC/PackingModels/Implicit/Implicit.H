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
    Foam::PackingModels::Implicit

Description
    Implicit model for applying an inter-particle stress to the particles.

    The time evolution of particulate volume fraction is solved for implicitly
    on the eulerian mesh. The computed flux is then applied to the lagrangian
    field. The gravity force can optionally be applied to the particles as part
    of this model, using the keyword "applyGravity".

SourceFiles
    Implicit.C

\*---------------------------------------------------------------------------*/

#ifndef Implicit_H
#define Implicit_H

#include "PackingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PackingModels
{

/*---------------------------------------------------------------------------*\
                          Class Implicit Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class Implicit
:
    public PackingModel<CloudType>
{
private:

    //- Private data
            
        //- Volume fraction field
        volScalarField alpha_;

        //- Correction flux
        tmp<surfaceScalarField> phiCorrect_;

        //- Correction cell-centred velocity
        tmp<volVectorField> uCorrect_;

        //- Flag to indicate whether gravity is applied
        Switch applyGravity_;

        //- Minimum stable volume fraction
        scalar alphaMin_;

        //- Minimum stable density
        scalar rhoMin_;


public:

    //- Runtime type information
    TypeName("implicit");

    // Constructors

        //- Construct from components
        Implicit(const dictionary& dict, CloudType& owner);

        //- Construct copy
        Implicit(const Implicit<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<PackingModel<CloudType> > clone() const
        {
            return autoPtr<PackingModel<CloudType> >
            (
                new Implicit<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~Implicit();


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
#   include "Implicit.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
