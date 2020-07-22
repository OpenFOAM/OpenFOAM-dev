/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2020 OpenFOAM Foundation
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
    Foam::DenseDragForce

\*---------------------------------------------------------------------------*/

#ifndef DenseDragForce_H
#define DenseDragForce_H

#include "ParticleForce.H"
#include "interpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class DenseDragForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class DenseDragForce
:
    public ParticleForce<CloudType>
{
protected:

    // Protected Data

        //- Name of the carrier volume fraction field
        const word alphacName_;

        //- The carrier volume fraction field
        autoPtr<volScalarField> alphacPtr_;

        //- Interpolation for the carrier volume fraction field
        autoPtr<interpolation<scalar>> alphacInterpPtr_;


public:

    // Constructors

        //- Construct from mesh
        DenseDragForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict,
            const word& typeName
        );

        //- Construct copy
        DenseDragForce(const DenseDragForce<CloudType>& df);


    //- Destructor
    virtual ~DenseDragForce();


    // Member Functions

        // Access

            //- Return the volume fraction interpolation
            const interpolation<scalar>& alphacInterp() const;

        // Evaluation

            //- Cache fields
            virtual void cacheFields(const bool store);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DenseDragForce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
