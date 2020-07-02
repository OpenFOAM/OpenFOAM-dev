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
    Foam::MPLICU

Description
    Velocity-weighted Multicut Piecewise-Linear Interface Calculation (MPLICU)
    corrected scheme is a surface interpolation scheme for flux calculation in
    advection of a bounded variable, e.g. phase fraction and for interface
    capturing in the volume of fluid (VoF) method.

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

    Additionally the face point velocity values are used to calculate the face
    flux which is likely to be more accurate in the presence of high shear.

    Example:
    \verbatim
    divSchemes
    {
    .
    .
    div(phi,alpha1)      Gauss MPLICU;
    .
    .
    }
    \endverbatim

See also
    Foam::MPLIC
    Foam::PLIC
    Foam::PLICU
    Foam::interfaceCompression

SourceFiles
    MPLICU.C

\*---------------------------------------------------------------------------*/

#ifndef MPLICU_H
#define MPLICU_H

#include "MPLIC.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class MPLICU Declaration
\*---------------------------------------------------------------------------*/

class MPLICU
:
    public MPLIC
{

public:

    //- Runtime type information
    TypeName("MPLICU");


    // Constructors

        //- Construct from faceFlux and Istream
        MPLICU
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            MPLIC(mesh, faceFlux, is)
        {}


    // Member Functions

        //- Return the interpolation weighting factors
        virtual tmp<surfaceScalarField> weights
        (
            const GeometricField<scalar, fvPatchField, volMesh>& vf
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
        void operator=(const MPLICU&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
