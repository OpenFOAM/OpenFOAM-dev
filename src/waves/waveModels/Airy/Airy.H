/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2021 OpenFOAM Foundation
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
    Foam::waveModels::Airy

Description
    First-order wave model.

    Reference:
    \verbatim
        Stokes, G.G. (1847)
        On the theory of oscillatory waves.
        Transactions of the Cambridge Philosophical Society, 8, 441.
    \endverbatim

    See the leading terms of equations 18 and 19.

SourceFiles
    Airy.C

\*---------------------------------------------------------------------------*/

#ifndef Airy_H
#define Airy_H

#include "waveModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveModels
{

/*---------------------------------------------------------------------------*\
                            Class Airy Declaration
\*---------------------------------------------------------------------------*/

class Airy
:
    public waveModel
{
    // Private Data

        //- Depth [m]
        const scalar depth_;

        //- Peak-to-peak mean amplitude [m]
        autoPtr<Function1<scalar>> amplitude_;

        //- Peak-to-peak length [m]
        const scalar length_;

        //- Phase offset [rad]
        const scalar phase_;


private:

    // Private Member Functions

        //- Read and return the wave length from the dictionary. Either reads
        //  the length directly, or reads the period and depth and calculates
        //  the length.
        static scalar length
        (
            const dictionary& dict,
            const scalar depth,
            const scalar amplitude,
            const scalar g,
            scalar (*celerityPtr)(scalar, scalar, scalar, scalar)
        );


protected:

    // Static Protected Member Functions

        //- The angular wavenumber [rad/m]
        static scalar k(const scalar length);

        //- Return whether shallow and intermediate effects are to be omitted
        static bool deep(const scalar length, const scalar depth);

        //- The wave celerity [m/s]
        static scalar celerity
        (
            const scalar depth,
            const scalar amplitude,
            const scalar length,
            const scalar g
        );


    // Protected Member Functions

        //- The angular wavenumber [rad/m]
        scalar k() const
        {
            return k(length_);
        }

        //- Return whether shallow and intermediate effects are to be omitted
        bool deep() const
        {
            return deep(depth(), length());
        }

        //- The wave celerity [m/s]
        virtual scalar celerity() const
        {
            return celerity(depth(), amplitude(), length(), g());
        }

        //- Angle of the oscillation [rad]
        tmp<scalarField> angle
        (
            const scalar t,
            const scalarField& x
        ) const;

        //- Return the non-dimensionalised i-th harmonic of the velocity
        tmp<vector2DField> vi
        (
            const label i,
            const scalar t,
            const vector2DField& xz
        ) const;


public:

    //- Runtime type information
    TypeName("Airy");


    // Constructors

        //- Construct a copy
        Airy(const Airy& wave);

        //- Construct from a dictionary and gravity
        Airy
        (
            const dictionary& dict,
            const scalar g,
            const word& modelName = Airy::typeName,
            scalar (*modelCelerity)(scalar, scalar, scalar, scalar) =
                &Airy::celerity
        );

        //- Construct a clone
        virtual autoPtr<waveModel> clone() const
        {
            return autoPtr<waveModel>(new Airy(*this));
        }


    //- Destructor
    virtual ~Airy();


    // Member Functions

        // Access

            //- Get the depth
            scalar depth() const
            {
                return depth_;
            }

            //- Get the amplitude
            scalar amplitude(const scalar t) const
            {
                return amplitude_->value(t);
            }

            //- Get the amplitude at steady state
            scalar amplitude() const
            {
                return amplitude_->value(great);
            }

            //- Get the length
            scalar length() const
            {
                return length_;
            }

            //- Get the phase
            scalar phase() const
            {
                return phase_;
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

        //- Write
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
