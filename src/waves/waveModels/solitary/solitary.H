/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2020 OpenFOAM Foundation
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
    Foam::waveModels::solitary

Description
    Solitary wave model.

    Reference:
    \verbatim
        Dean, R. G., & Dalrymple, R. A. (1991).
        Water wave mechanics for engineers and scientists (Vol. 2).
        World Scientific Publishing Company.
    \endverbatim

    See pages 314-317.

SourceFiles
    solitary.C

\*---------------------------------------------------------------------------*/

#ifndef solitary_H
#define solitary_H

#include "waveModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveModels
{

/*---------------------------------------------------------------------------*\
                          Class solitary Declaration
\*---------------------------------------------------------------------------*/

class solitary
:
    public waveModel
{
    // Private Data

        //- Offset [m]
        const scalar offset_;

        //- Depth [m]
        const scalar depth_;


    // Private Member Functions

        //- The wavenumber [1/m]
        scalar k(const scalar t) const;

        //- The dimensionless amplitude [1]
        scalar alpha(const scalar t) const;

        //- The wave celerity [m/s]
        scalar celerity(const scalar t) const;

        //- The evolution parameter [1]
        //  This is analogous to the oscillation angle of a periodic wave
        tmp<scalarField> parameter
        (
            const scalar t,
            const scalarField& x
        ) const;

        //- The dimensionless elevation [1]
        tmp<scalarField> Pi
        (
            const scalar t,
            const scalarField& x
        ) const;


public:

    //- Runtime type information
    TypeName("solitary");


    // Constructors

        //- Construct a copy
        solitary(const solitary& wave);

        //- Construct from a database and a dictionary
        solitary(const objectRegistry& db, const dictionary& dict);

        //- Construct a clone
        virtual autoPtr<waveModel> clone() const
        {
            return autoPtr<waveModel>(new solitary(*this));
        }


    //- Destructor
    virtual ~solitary();


    // Member Functions

        // Access

            //- Get the offset
            scalar offset() const
            {
                return offset_;
            }

            //- Get the depth
            scalar depth() const
            {
                return depth_;
            }

        //- Get the wave elevation at a given time and local coordinates. Local
        //  x is aligned with the direction of propagation.
        virtual tmp<scalarField> elevation
        (
            const scalar t,
            const scalarField& x
        ) const;

        //- Get the wave velocity at a given time and local coordinates. Local
        //  x is aligned with the direction of propagation, and z with negative
        //  gravity.
        virtual tmp<vector2DField> velocity
        (
            const scalar t,
            const vector2DField& xz
        ) const;

        //- Get the wave pressure at a given time and local coordinates. Local
        //  x is aligned with the direction of propagation, and z with negative
        //  gravity.
        virtual tmp<scalarField> pressure
        (
            const scalar t,
            const vector2DField& xz
        ) const;

        //- Write
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
