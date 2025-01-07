/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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
    Foam::CarrierField

Description
    A field interpolated from the carrier to the cloud. Uses CloudDerivedField
    to provide flexible access and caching.

See also
    Foam::CloudDerivedField

SourceFiles
    CarrierField.C

\*---------------------------------------------------------------------------*/

#ifndef CarrierField_H
#define CarrierField_H

#include "CloudDerivedField.H"
#include "interpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class CarrierField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class CarrierField
:
    public CloudDerivedField<Type>
{
    // Private Classes

        //- Class to virtualise an evaluation procedure
        class Functor;

        //- Class to store an evaluation function
        template<class F>
        class Function;


    // Private Data

        //- The carrier field name
        const word name_;

        //- The function
        const autoPtr<Functor> functorPtr_;

        //- Field reference or cache
        mutable tmp<VolField<Type>> tpsi_;

        //- Interpolation engine
        mutable autoPtr<interpolation<Type>> psiInterpolationPtr_;

        //- Old-time interpolation engine
        mutable autoPtr<interpolation<Type>> psi0InterpolationPtr_;


    // Private Member Functions

        //- Interpolate and return the result
        tmp<LagrangianSubField<Type>> interpolate
        (
            const LagrangianModelRef&,
            const LagrangianSubMesh& subMesh
        ) const;


public:

    // Constructors

        //- Construct from a reference to a carrier field
        CarrierField(const VolField<Type>&);

        //- Construct from a name and an evaluation function
        template<class F>
        CarrierField(const word& name, const F& f);


    // Member Functions

        //- Access the carrier field
        const VolField<Type>& psi() const;

        //- Carrier field name
        const word& name() const;

        //- Reset
        void reset(const bool predict);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CarrierField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
