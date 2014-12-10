/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "coordinateRotation.H"
#include "dictionary.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(coordinateRotation, 0);
    defineRunTimeSelectionTable(coordinateRotation, dictionary);
    defineRunTimeSelectionTable(coordinateRotation, objectRegistry);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::symmTensor Foam::coordinateRotation::transformPrincipal
(
    const tensor& tt,
    const vector& st
) const
{
    return symmTensor
    (
        tt.xx()*st.x()*tt.xx()
      + tt.xy()*st.y()*tt.xy()
      + tt.xz()*st.z()*tt.xz(),

        tt.xx()*st.x()*tt.yx()
      + tt.xy()*st.y()*tt.yy()
      + tt.xz()*st.z()*tt.yz(),

        tt.xx()*st.x()*tt.zx()
      + tt.xy()*st.y()*tt.zy()
      + tt.xz()*st.z()*tt.zz(),

        tt.yx()*st.x()*tt.yx()
      + tt.yy()*st.y()*tt.yy()
      + tt.yz()*st.z()*tt.yz(),

        tt.yx()*st.x()*tt.zx()
      + tt.yy()*st.y()*tt.zy()
      + tt.yz()*st.z()*tt.zz(),

        tt.zx()*st.x()*tt.zx()
      + tt.zy()*st.y()*tt.zy()
      + tt.zz()*st.z()*tt.zz()
    );

}

// ************************************************************************* //
