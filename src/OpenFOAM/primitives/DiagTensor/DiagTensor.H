/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::DiagTensor

Description
    Templated 3D DiagTensor derived from VectorSpace.

    Adding construction from 3 components, element access using xx(), yy()
    and zz() member functions and the inner-product (dot-product) and
    outer-product operators.

SourceFiles
    DiagTensorI.H

\*---------------------------------------------------------------------------*/

#ifndef DiagTensor_H
#define DiagTensor_H

#include "Tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class DiagTensor Declaration
\*---------------------------------------------------------------------------*/

template<class Cmpt>
class DiagTensor
:
    public VectorSpace<DiagTensor<Cmpt>, Cmpt, 3>
{

public:

    //- Equivalent type of labels used for valid component indexing
    typedef DiagTensor<label> labelType;


    // Member constants

        //- Rank of DiagTensor is 2
        static const direction rank = 2;


    //- Component labeling enumeration
    enum components { XX, YY, ZZ };


    // Constructors

        //- Construct null
        inline DiagTensor();

        //- Construct initialized to zero
        inline DiagTensor(const Foam::zero);

        //- Construct given VectorSpace
        template<class Cmpt2>
        inline DiagTensor(const VectorSpace<DiagTensor<Cmpt2>, Cmpt2, 3>&);

        //- Construct given three components
        inline DiagTensor(const Cmpt& txx, const Cmpt& tyy, const Cmpt& tzz);

        //- Construct from Istream
        inline DiagTensor(Istream&);


    // Member Functions

        // Access

            inline const Cmpt& xx() const;
            inline const Cmpt& yy() const;
            inline const Cmpt& zz() const;

            inline Cmpt& xx();
            inline Cmpt& yy();
            inline Cmpt& zz();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Include inline implementations
#include "DiagTensorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
