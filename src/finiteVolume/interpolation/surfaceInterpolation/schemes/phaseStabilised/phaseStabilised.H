/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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
    Foam::phaseStabilised

Description
    Phase-stabilised interpolation scheme

    which applies upwind to the faces at which the upwind phase-fraction is less
    than 1e-3.

SourceFiles
    phaseStabilised.C

\*---------------------------------------------------------------------------*/

#ifndef phaseStabilised_H
#define phaseStabilised_H

#include "upwind.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class phaseStabilised Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class phaseStabilised
:
    public surfaceInterpolationScheme<Type>
{
    // Private member data

        const surfaceScalarField& faceFlux_;
        tmp<surfaceInterpolationScheme<Type>> tScheme_;

        upwind<scalar> upwind_;

        const volScalarField& alpha_;


    // Private member functions

        const volScalarField& lookupAlpha
        (
            const surfaceScalarField& faceFlux
        ) const;

        tmp<surfaceScalarField> lambdaf() const;


public:

    //- Runtime type information
    TypeName("phaseStabilised");


    // Constructors

        //- Construct from mesh and Istream
        phaseStabilised(const fvMesh& mesh, Istream& is);


        //- Construct from mesh, faceFlux and Istream
        phaseStabilised
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        );

        //- Disallow default bitwise copy construction
        phaseStabilised(const phaseStabilised&) = delete;


    // Member Functions

        //- Return the interpolation weighting factors
        tmp<surfaceScalarField> weights
        (
            const VolField<Type>& vf
        ) const;

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return tScheme_().corrected();
        }

        //- Return the explicit correction to the face-interpolate
        //  set to zero on the near-boundary faces where upwind is applied
        virtual tmp<SurfaceField<Type>>
        correction
        (
            const VolField<Type>& vf
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const phaseStabilised&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "phaseStabilisedTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
