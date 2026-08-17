/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2026 OpenFOAM Foundation
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
    Foam::fv::filmGaussGrad

Description
    Second-order Gauss gradient scheme for film

SourceFiles
    filmGaussGrad.C

\*---------------------------------------------------------------------------*/

#ifndef filmGaussGrad_H
#define filmGaussGrad_H

#include "gaussGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                      Class filmGaussGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class filmGaussGrad
:
    public gaussGrad<Type>
{
public:

    //- Runtime type information
    TypeName("filmGauss");


    // Constructors

        //- Construct from mesh
        filmGaussGrad(const fvMesh& mesh)
        :
            gaussGrad<Type>(mesh)
        {}

        //- Construct from Istream
        filmGaussGrad(const fvMesh& mesh, Istream& is)
        :
            gaussGrad<Type>(mesh, is)
        {}

        //- Disallow default bitwise copy construction
        filmGaussGrad(const filmGaussGrad&) = delete;


    // Member Functions

        //- Set the grad field to the gradient of the volField vf
        virtual void calcGrad
        (
            VolInternalField<typename outerProduct<vector, Type>::type>& grad,
            const VolField<Type>& vf
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const filmGaussGrad&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "filmGaussGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
