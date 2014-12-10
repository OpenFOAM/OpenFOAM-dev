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
    Foam::aspectRatioModels::Wellek

Description
    Aspect ratio model of Wellek et al.

    References:
    \verbatim
        "Implementation and Comparison of Correlations for interfacial Forces
        in a Gas-Liquid System within an Euler-Euler Framework"
        M Otromke
        PhD Thesis
        April 2013
    \endverbatim

    \verbatim
        "Shape of liquid drops moving in liquid media"
        R M Wellek, A K Agrawal, A H P Skelland
        International Journal of Multiphase Flow
        Volume 12, Issue 5, Pages 854-862, September 1966
    \endverbatim

SourceFiles
    Wellek.C

\*---------------------------------------------------------------------------*/

#ifndef Wellek_H
#define Wellek_H

#include "aspectRatioModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace aspectRatioModels
{

/*---------------------------------------------------------------------------*\
                       Class Wellek Declaration
\*---------------------------------------------------------------------------*/

class Wellek
:
    public aspectRatioModel
{
public:

    //- Runtime type information
    TypeName("Wellek");


    // Constructors

        //- Construct from a dictionary and an ordered phase pair
        Wellek
        (
            const dictionary& dict,
            const orderedPhasePair& pair
        );


    //- Destructor
    virtual ~Wellek();


    // Member Functions

        //- Aspect ratio
        virtual tmp<volScalarField> E() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace aspectRatioModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
