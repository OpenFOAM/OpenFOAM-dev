/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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
    Foam::MPLIC

Description
    Multicut Piecewise-Linear Interface Calculation (MPLIC) corrected scheme is
    a surface interpolation scheme for flux calculation in advection of a
    bounded variable, e.g. phase fraction and for interface capturing in the
    volume of fluid (VoF) method.

    The interface is represented by multiple cuts which split each cell to match
    the volume fraction of the phase in the cell. The cut planes are oriented
    according to the point field of the local phase fraction.  The phase
    fraction at each cell face - the interpolated value - is then calculated
    from the face area on either side of the cuts.

    Three progressively more complex algorithms are used to ensure the cell
    volume fraction is accurately reproduced:
      -# single cut: cuts all the cell faces regardless the order
      -# multi cut: topological face-edge-face walk which can split cell into
         multiple sub-volumes
      -# tetrahedron cut: decomposes cell into tetrahedrons which are cut

    Example:
    \verbatim
    divSchemes
    {
        .
        .
        div(phi,alpha)      Gauss MPLIC;
        .
        .
    }
    \endverbatim

See also
    Foam::MPLICU
    Foam::PLIC
    Foam::PLICU
    Foam::interfaceCompression

SourceFiles
    MPLIC.C

\*---------------------------------------------------------------------------*/

#ifndef MPLIC_H
#define MPLIC_H

#include "surfaceInterpolationScheme.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class MPLIC Declaration
\*---------------------------------------------------------------------------*/

class MPLIC
:
    public surfaceInterpolationScheme<scalar>
{
protected:

    // Private member data

        const surfaceScalarField& phi_;


    // Protected member functions

        //- Set alphaPhi for the faces of the given cell
        void setCellAlphaf
        (
            const label celli,
            const scalarField& phi,
            scalarField& alphaf,
            boolList& correctedFaces,
            const DynamicList<scalar>& cellAlphaf,
            const fvMesh& mesh
        ) const;

        //- Return alpha interpolation
        tmp<surfaceScalarField> surfaceAlpha
        (
            const volScalarField& alpha,
            const surfaceScalarField& phi,
            scalarField& spicedTvff,
            const bool unweighted,
            const scalar tol,
            const bool isMPLIC = true
        ) const;


public:

    //- Runtime type information
    TypeName("MPLIC");


    // Constructors

        //- Construct from faceFlux and Istream
        MPLIC
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            surfaceInterpolationScheme<scalar>(mesh),
            phi_(faceFlux)
        {}


    // Member Functions

        //- Return the interpolation weighting factors
        virtual tmp<surfaceScalarField> weights
        (
            const GeometricField<scalar, fvPatchField, volMesh>&
        ) const
        {
            NotImplemented;

            return tmp<surfaceScalarField>(nullptr);
        }

        //- Return the face-interpolate of the given cell field
        virtual tmp<surfaceScalarField> interpolate
        (
            const GeometricField<scalar, fvPatchField, volMesh>& vf
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const MPLIC&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
