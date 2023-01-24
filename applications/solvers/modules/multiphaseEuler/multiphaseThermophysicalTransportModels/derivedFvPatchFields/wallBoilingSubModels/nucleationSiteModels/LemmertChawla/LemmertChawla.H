/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2023 OpenFOAM Foundation
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
    Foam::wallBoilingModels::nucleationSiteModels::LemmertChawla

Description
    Lemmert & Chawla function for nucleation site density,
    correlation by Egorov & Menter.

    References:
    \verbatim
        Lemmert, M., & Chawla, J. M. (1977).
        Influence of flow velocity on surface boiling heat transfer coefficient.
        Heat Transfer in Boiling, 237, 247.

        Egorov, Y., & Menter, F. (2004).
        Experimental implementation of the RPI wall boiling model in CFX-5.6.
        Staudenfeldweg, 12, 83624.
        Technical Report ANSYS/TR-04-10, ANSYS Gmbh.
    \endverbatim

SourceFiles
    LemmertChawla.C

\*---------------------------------------------------------------------------*/

#ifndef LemmertChawla_H
#define LemmertChawla_H

#include "nucleationSiteModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace nucleationSiteModels
{

/*---------------------------------------------------------------------------*\
                        Class LemmertChawla Declaration
\*---------------------------------------------------------------------------*/

class LemmertChawla
:
    public nucleationSiteModel
{
    // Private Data

        //- Coefficient for nucleation site density
        scalar Cn_;

        //- Reference nucleations density
        scalar NRef_;

        //- Reference temperature difference
        scalar deltaTRef_;


public:

    //- Runtime type information
    TypeName("LemmertChawla");

    // Constructors

        //- Construct from a dictionary
        LemmertChawla(const dictionary& dict);

        //- Copy construct
        LemmertChawla(const LemmertChawla& model);

        //- Construct and return a clone
        virtual autoPtr<nucleationSiteModel> clone() const
        {
            return autoPtr<nucleationSiteModel>(new LemmertChawla(*this));
        }


    //- Destructor
    virtual ~LemmertChawla();


    // Member Functions

        //- Calculate and return the nucleation-site density
        virtual tmp<scalarField> nucleationSiteDensity
        (
            const phaseModel& liquid,
            const phaseModel& vapour,
            const label patchi,
            const scalarField& Tl,
            const scalarField& Tsatw,
            const scalarField& L,
            const scalarField& dDep,
            const scalarField& fDep
        ) const;

        //- Write
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace nucleationSiteModels
} // End namespace wallBoilingModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
