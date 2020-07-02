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
    Foam::PLICU

Description
    Velocity-weighted Piecewise-Linear Interface Calculation (PLICU) corrected
    scheme is a surface interpolation scheme for flux calculation in advection
    of a bounded variable, e.g. phase fraction and for interface capturing in
    the volume of fluid (VoF) method.

    The interface is represented by single cuts which split each cell to match
    the volume fraction of the phase in the cell. The cut planes are oriented
    according to the point field of the local phase fraction.  The phase
    fraction at each cell face - the interpolated value - is then calculated
    from the face area on either side of the cut.  For cases where the
    single-cut does not accurately represent the cell volume fraction the
    specified default scheme is used, e.g. interfaceCompression.

    Additionally the face point velocity values are used to calculate the face
    flux which is likely to be more accurate in the presence of high shear.

    Example:
    \verbatim
    divSchemes
    {
        .
        .
        div(phi,alpha)     Gauss PLICU interfaceCompression vanLeer 1;
        .
        .
    }
    \endverbatim

See also
    Foam::PLIC
    Foam::MPLIC
    Foam::MPLICU
    Foam::interfaceCompression

SourceFiles
    PLICU.C

\*---------------------------------------------------------------------------*/

#ifndef PLICU_H
#define PLICU_H

#include "PLIC.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class PLIC Declaration
\*---------------------------------------------------------------------------*/

class PLICU
:
    public PLIC
{

public:

    //- Runtime type information
    TypeName("PLICU");


    // Constructors

        //- Construct from faceFlux and Istream
        PLICU
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            PLIC(mesh, faceFlux, is)
        {}


    // Member Functions

        //- Return the face-interpolate of the given cell field
        virtual tmp<surfaceScalarField> interpolate
        (
            const volScalarField& vf
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const PLICU&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
