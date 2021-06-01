/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2021 OpenFOAM Foundation
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
    Foam::Barycentric

Description
    Templated 3D Barycentric derived from VectorSpace. Has 4 components, one of
    which is redundant.

SourceFiles
    BarycentricI.H

\*---------------------------------------------------------------------------*/

#ifndef Barycentric_H
#define Barycentric_H

#include "VectorSpace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class Barycentric Declaration
\*---------------------------------------------------------------------------*/

template<class Cmpt>
class Barycentric
:
    public VectorSpace<Barycentric<Cmpt>, Cmpt, 4>
{
public:

    //- Equivalent type of labels used for valid component indexing
    typedef Barycentric<label> labelType;


    // Member constants

        //- Rank of Barycentric is 1
        static const direction rank = 1;


    //- Component labeling enumeration
    enum components { A, B, C, D };


    // Constructors

        //- Construct null
        inline Barycentric();

        //- Construct initialised to zero
        inline Barycentric(const Foam::zero);

        //- Construct given four components
        inline Barycentric
        (
            const Cmpt& va,
            const Cmpt& vb,
            const Cmpt& vc,
            const Cmpt& vd
        );


    // Member Functions

        // Access

            inline const Cmpt& a() const;
            inline const Cmpt& b() const;
            inline const Cmpt& c() const;
            inline const Cmpt& d() const;

            inline Cmpt& a();
            inline Cmpt& b();
            inline Cmpt& c();
            inline Cmpt& d();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "BarycentricI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
