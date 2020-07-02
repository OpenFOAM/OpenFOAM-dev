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
    Foam::PLIC

Description
    Piecewise-Linear Interface Calculation (PLIC) corrected scheme is a surface
    interpolation scheme for flux calculation in advection of a bounded
    variable, e.g. phase fraction and for interface capturing in the volume of
    fluid (VoF) method.

    The interface is represented by single cuts which split each cell to match
    the volume fraction of the phase in the cell. The cut planes are oriented
    according to the point field of the local phase fraction.  The phase
    fraction at each cell face - the interpolated value - is then calculated
    from the face area on either side of the cut.  For cases where the
    single-cut does not accurately represent the cell volume fraction the
    specified default scheme is used, e.g. interfaceCompression.

    Example:
    \verbatim
    divSchemes
    {
        .
        .
        div(phi,alpha)     Gauss PLIC interfaceCompression vanLeer 1;
        .
        .
    }
    \endverbatim

See also
    Foam::PLICU
    Foam::MPLIC
    Foam::MPLICU
    Foam::interfaceCompression

SourceFiles
    PLIC.C

\*---------------------------------------------------------------------------*/

#ifndef PLIC_H
#define PLIC_H

#include "MPLIC.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class PLIC Declaration
\*---------------------------------------------------------------------------*/

class PLIC
:
    public MPLIC
{

protected:

    // Protected member data

        const surfaceScalarField& phi_;

        //- Base scheme to which the compression is applied
        tmp<surfaceInterpolationScheme<scalar>> tScheme_;


public:

    //- Runtime type information
    TypeName("PLIC");


    // Constructors

        //- Construct from faceFlux and Istream
        PLIC
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            MPLIC(mesh, faceFlux, is),
            phi_(faceFlux),
            tScheme_
            (
                surfaceInterpolationScheme<scalar>::New(mesh, faceFlux, is)
            )
        {}


    // Member Functions

        //- Return the face-interpolate of the given cell field
        virtual tmp<surfaceScalarField> interpolate
        (
            const volScalarField& vf
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const PLIC&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
