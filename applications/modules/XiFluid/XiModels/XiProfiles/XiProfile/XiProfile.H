/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024-2025 OpenFOAM Foundation
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
    Foam::XiProfile

Description
    Base class for flame wrinkling profiles

    Flame wrinkling profiles are functions of \c b only for simplicity and
    numerical stability.  The profiles are usually formulated so that the flame
    wrinkling \c Xi is equal to the equilibrium value in the middle of the
    flame, i.e. where \c b = 0.5, but this does not guarantee that the flame
    propagates at the correct speed with heat release due to the slight
    asymmetry in \c b distribution caused by the heat release dilatation.

    In principle a better distribution function for \c Xi would be based
    directly on the \c laplacian(b) such that \c Xi equals the equilibrium
    value at the point in the flame where the Laplacian of \c b is 0, e.g.
    \verbatim
        Xi = max(XiEq + 0.1*XiEqrho*fvc::laplacian(Db, b)/(rhou*Db*sqr(mgb)), 1)
    \endverbatim
    While this does provide the correct speed for a well resolved 1D flame it is
    prone to cause serious numerical instability and is not applicable to real
    flame simulations.

    A more practical approach to correct the flame speed for the dilatation
    effect would be to introduce a small amount of off-centring to the \c Xi
    distribution as a function of \c b which would need to be a function of at
    least the density ratio and possibly the equilibrium \c Xi also -- this is
    the subject of current research.

SourceFiles
    XiProfile.C

\*---------------------------------------------------------------------------*/

#ifndef XiProfile_H
#define XiProfile_H

#include "volFields.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class XiProfile Declaration
\*---------------------------------------------------------------------------*/

class XiProfile
{

protected:

    // Protected data

        const volScalarField& b_;


    // Protected member functions

        //- Update coefficients from given dictionary
        virtual bool readCoeffs(const dictionary& dict) = 0;


public:

    //- Runtime type information
    TypeName("XiProfile");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            XiProfile,
            dictionary,
            (
                const dictionary& dict,
                const volScalarField& b
            ),
            (
                dict,
                b
            )
        );


    // Constructors

        //- Construct from components
        XiProfile(const volScalarField& b);

        //- Disallow default bitwise copy construction
        XiProfile(const XiProfile&) = delete;


    // Selectors

        //- Return a reference to the selected XiEq model
        static autoPtr<XiProfile> New
        (
            const dictionary& XiProperties,
            const volScalarField& b
        );


    //- Destructor
    virtual ~XiProfile();


    // Member Functions

        //- Return the flame-wrinkling profile
        virtual tmp<volScalarField> profile() const = 0;

        //- Update properties from the given XiProperties dictionary
        bool read(const dictionary& XiProperties);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const XiProfile&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
