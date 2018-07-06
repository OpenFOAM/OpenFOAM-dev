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
    Foam::PhiScheme

Description
    Class to create the weighting-factors based on the face-flux.

    The particular differencing scheme class is supplied as a template
    argument, the weight function of which is called by the weight function
    of this class for the internal faces as well as faces of coupled
    patches (e.g. processor-processor patches). The weight function is
    supplied with the central-differencing weighting factor, the face-flux,
    the face neighbour cell values and the face area.

    This code organisation is both neat and efficient, allowing for
    convenient implementation of new schemes to run on parallelised cases.

SourceFiles
    PhiScheme.C

\*---------------------------------------------------------------------------*/

#ifndef PhiScheme_H
#define PhiScheme_H

#include "limitedSurfaceInterpolationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class PhiScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class PhiLimiter>
class PhiScheme
:
    public limitedSurfaceInterpolationScheme<Type>,
    public PhiLimiter
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        PhiScheme(const PhiScheme&);

        //- Disallow default bitwise assignment
        void operator=(const PhiScheme&);


public:

    //- Runtime type information
    TypeName("PhiScheme");


    // Constructors

        //- Construct from mesh, faceFlux and blendingFactor
        PhiScheme
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            const PhiLimiter& weight
        )
        :
            limitedSurfaceInterpolationScheme<Type>(mesh, faceFlux),
            PhiLimiter(weight)
        {}

        //- Construct from mesh and Istream.
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        PhiScheme
        (
            const fvMesh& mesh,
            Istream& is
        )
        :
            limitedSurfaceInterpolationScheme<Type>(mesh, is),
            PhiLimiter(is)
        {}

        //- Construct from mesh, faceFlux and Istream
        PhiScheme
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            limitedSurfaceInterpolationScheme<Type>(mesh, faceFlux),
            PhiLimiter(is)
        {}


    // Member Functions

        //- Return the interpolation weighting factors
        virtual tmp<surfaceScalarField> limiter
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makePhiSurfaceInterpolationScheme(SS, WEIGHT, TYPE)                    \
                                                                               \
typedef PhiScheme<TYPE, WEIGHT> Phischeme##WEIGHT_;                            \
defineTemplateTypeNameAndDebugWithName(Phischeme##WEIGHT_, #SS, 0);            \
                                                                               \
surfaceInterpolationScheme<TYPE>::addMeshConstructorToTable                    \
<PhiScheme<TYPE, WEIGHT>> add##SS##TYPE##MeshConstructorToTable_;              \
                                                                               \
surfaceInterpolationScheme<TYPE>::addMeshFluxConstructorToTable                \
<PhiScheme<TYPE, WEIGHT>> add##SS##TYPE##MeshFluxConstructorToTable_;          \
                                                                               \
limitedSurfaceInterpolationScheme<TYPE>::addMeshConstructorToTable             \
<PhiScheme<TYPE, WEIGHT>> add##SS##TYPE##MeshConstructorToLimitedTable_;       \
                                                                               \
limitedSurfaceInterpolationScheme<TYPE>::addMeshFluxConstructorToTable         \
<PhiScheme<TYPE, WEIGHT>> add##SS##TYPE##MeshFluxConstructorToLimitedTable_;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PhiScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
