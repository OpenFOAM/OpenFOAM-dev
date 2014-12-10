/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 OpenFOAM Foundation
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
    Foam::dragModels::Lain

Description
    Drag model of Lain et al.

    References:
    \verbatim
        "Implementation and Comparison of Correlations for interfacial Forces
        in a Gas-Liquid System within an Euler-Euler Framework"
        M Otromke
        PhD Thesis
        April 2013
    \endverbatim

    \verbatim
        "Modelling hydrodynamics and turbulence in a bubble column using the
        Euler-Lagrange procedure"
        S Lain, D Brodera, M Sommerfelda, M F Goza
        International Journal of Multiphase Flow
        Volume 28, Issue 8, Pages 1381-1407, August 2002
    \endverbatim

SourceFiles
    Lain.C

\*---------------------------------------------------------------------------*/

#ifndef Lain_H
#define Lain_H

#include "dragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace dragModels
{

/*---------------------------------------------------------------------------*\
                           Class Lain Declaration
\*---------------------------------------------------------------------------*/

class Lain
:
    public dragModel
{
public:

    //- Runtime type information
    TypeName("Lain");


    // Constructors

        //- Construct from a dictionary and a phase pair
        Lain
        (
            const dictionary& dict,
            const phasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~Lain();


    // Member Functions

        //- Drag coefficient
        virtual tmp<volScalarField> CdRe() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dragModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
