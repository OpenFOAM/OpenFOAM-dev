/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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
    Foam::thirdBodyEfficiencies

Description
    Third body efficiencies

SourceFiles
    thirdBodyEfficienciesI.H

\*---------------------------------------------------------------------------*/

#ifndef thirdBodyEfficiencies_H
#define thirdBodyEfficiencies_H

#include "scalarField.H"
#include "speciesTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class thirdBodyEfficiencies;
Ostream& operator<<(Ostream&, const thirdBodyEfficiencies&);


/*---------------------------------------------------------------------------*\
                    Class thirdBodyEfficiencies Declaration
\*---------------------------------------------------------------------------*/

class thirdBodyEfficiencies
{
    // Private Data

        const speciesTable& species_;

        //- Default efficiency
        scalar defaultEfficiency_;

        //- Third-body efficiencies
        scalarField efficiencies_;

        //- Specie indices with non-unity third-body efficiencies
        labelList thirdBodySpecies_;

        //- Third-body efficiencies - reference
        scalarField efficienciesMdefault_;


public:

    // Constructors

        //- Construct from table of efficiencies
        inline thirdBodyEfficiencies
        (
            const speciesTable& species,
            const scalar defaultEfficiency,
            const List<Tuple2<word, scalar>>& efficiencies
        );

        //- Construct from dictionary
        inline thirdBodyEfficiencies
        (
            const speciesTable& species,
            const dictionary& dict
        );


    // Member Functions

        //- Calculate and return M, the concentration of the third-bodies
        inline scalar M(const reaction::state& s) const;

        //- Return the derivative of M, w.r.t. the species concentrations
        inline const scalarField& dMdc() const;

        //- Write to stream
        inline void write(Ostream& os) const;



    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const thirdBodyEfficiencies&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "thirdBodyEfficienciesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
