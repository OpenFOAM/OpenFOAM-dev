/*---------------------------------------------------------------------------*\
 =========                   |
 \\      /   F ield          | OpenFOAM: The Open Source CFD Toolbox
  \\    /    O peration      |
   \\  /     A nd            | Copyright (C) 2012-2018 OpenFOAM Foundation
    \\/      M anipulation   |
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

Application
    foamyHexMeshSurfaceSimplify

Description
    Simplifies surfaces by resampling.

    Uses Thomas Lewiner's topology preserving MarchingCubes.
    (http://zeus.mat.puc-rio.br/tomlew/tomlew_uk.php)

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "searchableSurfaces.H"
#include "conformationSurfaces.H"
#include "triSurfaceMesh.H"

#include "opt_octree.h"
#include "cube.h"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class pointConversion
{
    const vector scale_;

    const vector offset_;

public:

    //- Construct from components
    pointConversion
    (
        const vector scale,
        const vector offset
    )
    :
        scale_(scale),
        offset_(offset)
    {}

    inline Point toLocal(const Foam::point& pt) const
    {
        Foam::point p = cmptMultiply(scale_, (pt + offset_));
        return Point(p.x(), p.y(), p.z());
    }

    inline Foam::point toGlobal(const Point& pt) const
    {
        point p(pt.x(), pt.y(), pt.z());
        return cmptDivide(p, scale_) - offset_;
    }
};




// For use in Fast-Dual Octree from Thomas Lewiner
class distanceCalc
:
    public ::data_access
{

    const Level min_level_;

    const conformationSurfaces& geometryToConformTo_;

    const pointConversion& converter_;


    // Private Member Functions

    scalar signedDistance(const Foam::point& pt) const
    {
        const searchableSurfaces& geometry = geometryToConformTo_.geometry();
        const labelList& surfaces = geometryToConformTo_.surfaces();

        static labelList nearestSurfaces;
        static scalarField distance;

        static pointField samples(1);
        samples[0] = pt;

        searchableSurfacesQueries::signedDistance
        (
            geometry,
            surfaces,
            samples,
            scalarField(1, great),
            volumeType::OUTSIDE,
            nearestSurfaces,
            distance
        );

        return distance[0];
    }


public:

    // Constructors

        //- Construct from components
        distanceCalc
        (
            Level max_level_,
            real iso_val_,
            Level min_level,
            const conformationSurfaces& geometryToConformTo,
            const pointConversion& converter
        )
        :
            data_access(max_level_,iso_val_),
            min_level_(min_level),
            geometryToConformTo_(geometryToConformTo),
            converter_(converter)
        {}


    //- Destructor
    virtual ~distanceCalc()
    {}


    // Member Functions

        //- Test function
        virtual bool need_refine( const Cube &c )
        {
            int l = c.lv() ;

            if ( l >= _max_level ) return false;
            if ( l < min_level_ ) return true;

            treeBoundBox bb
            (
                converter_.toGlobal
                (
                    Point
                    (
                        c.xmin(),
                        c.ymin(),
                        c.zmin()
                    )
                ),
                converter_.toGlobal
                (
                    Point
                    (
                        c.xmax(),
                        c.ymax(),
                        c.zmax()
                    )
                )
            );

            const searchableSurfaces& geometry =
                geometryToConformTo_.geometry();
            const labelList& surfaces =
                geometryToConformTo_.surfaces();


            //- Uniform refinement around surface
            {
                forAll(surfaces, i)
                {
                    if (geometry[surfaces[i]].overlaps(bb))
                    {
                        return true;
                    }
                }
                return false;
            }


            ////- Surface intersects bb (but not using intersection test)
            // scalar ccDist = signedDistance(bb.midpoint());
            // scalar ccVal = ccDist - _iso_val;
            // if (mag(ccVal) < small)
            //{
            //    return true;
            //}
            // const pointField points(bb.points());
            // forAll(points, pointi)
            //{
            //    scalar pointVal = signedDistance(points[pointi]) - _iso_val;
            //    if (ccVal*pointVal < 0)
            //    {
            //        return true;
            //    }
            //}
            // return false;


            ////- Refinement based on intersection with multiple planes.
            ////  Does not work well - too high a ratio between
            ////  neighbouring cubes.
            // const pointField points(bb.points());
            // const edgeList& edges = treeBoundBox::edges;
            // pointField start(edges.size());
            // pointField end(edges.size());
            // forAll(edges, i)
            //{
            //    start[i] = points[edges[i][0]];
            //    end[i] = points[edges[i][1]];
            //}
            // Foam::List<Foam::List<pointIndexHit>> hitInfo;
            // labelListList hitSurfaces;
            // searchableSurfacesQueries::findAllIntersections
            //(
            //    geometry,
            //    surfaces,
            //    start,
            //    end,
            //    hitSurfaces,
            //    hitInfo
            //);
            //
            //// Count number of intersections
            // label nInt = 0;
            // forAll(hitSurfaces, edgeI)
            //{
            //    nInt += hitSurfaces[edgeI].size();
            //}
            //
            // if (nInt == 0)
            //{
            //    // No surface intersected. See if there is one inside
            //    forAll(surfaces, i)
            //    {
            //        if (geometry[surfaces[i]].overlaps(bb))
            //        {
            //            return true;
            //        }
            //    }
            //    return false;
            //}
            //
            //// Check multiple surfaces
            // label baseSurfI = -1;
            // forAll(hitSurfaces, edgeI)
            //{
            //    const labelList& hSurfs = hitSurfaces[edgeI];
            //    forAll(hSurfs, i)
            //    {
            //        if (baseSurfI == -1)
            //        {
            //            baseSurfI = hSurfs[i];
            //        }
            //        else if (baseSurfI != hSurfs[i])
            //        {
            //            // Multiple surfaces
            //            return true;
            //        }
            //    }
            //}
            //
            //// Get normals
            // DynamicList<pointIndexHit> baseInfo(nInt);
            // forAll(hitInfo, edgeI)
            //{
            //    const Foam::List<pointIndexHit>& hits = hitInfo[edgeI];
            //    forAll(hits, i)
            //    {
            //        (void)hits[i].hitPoint();
            //        baseInfo.append(hits[i]);
            //    }
            //}
            // vectorField normals;
            // geometry[surfaces[baseSurfI]].getNormal(baseInfo, normals);
            // for (label i = 1; i < normals.size(); ++i)
            //{
            //    if ((normals[0] & normals[i]) < 0.9)
            //    {
            //        return true;
            //    }
            //}
            // labelList regions;
            // geometry[surfaces[baseSurfI]].getRegion(baseInfo, regions);
            // for (label i = 1; i < regions.size(); ++i)
            //{
            //    if (regions[0] != regions[i])
            //    {
            //        return true;
            //    }
            //}
            // return false;



            // samples[0] = point(c.xmin(), c.ymin(), c.zmin());
            // samples[1] = point(c.xmax(), c.ymin(), c.zmin());
            // samples[2] = point(c.xmax(), c.ymax(), c.zmin());
            // samples[3] = point(c.xmin(), c.ymax(), c.zmin());
            //
            // samples[4] = point(c.xmin(), c.ymin(), c.zmax());
            // samples[5] = point(c.xmax(), c.ymin(), c.zmax());
            // samples[6] = point(c.xmax(), c.ymax(), c.zmax());
            // samples[7] = point(c.xmin(), c.ymax(), c.zmax());

            // scalarField nearestDistSqr(8, great);
            //
            // Foam::List<pointIndexHit> nearestInfo;
            // surf_.findNearest(samples, nearestDistSqr, nearestInfo);
            // vectorField normals;
            // surf_.getNormal(nearestInfo, normals);
            //
            // for (label i = 1; i < normals.size(); ++i)
            //{
            //    if ((normals[0] & normals[i]) < 0.5)
            //    {
            //        return true;
            //    }
            //}
            // return false;

            //// Check if surface octree same level
            // const labelList elems(surf_.tree().findBox(bb));
            //
            // if (elems.size() > 1)
            //{
            //    return true;
            //}
            // else
            //{
            //  return false;
            //}
        }

        //- Data function
        virtual real value_at( const Cube &c )
        {
            return signedDistance(converter_.toGlobal(c)) - _iso_val;
        }
};


// Main program:

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Re-sample surfaces used in foamyHexMesh operation"
    );
    argList::validArgs.append("outputName");

    #include "setRootCase.H"
    #include "createTime.H"
    runTime.functionObjects().off();

    const fileName exportName = args.args()[1];

    Info<< "Reading surfaces as specified in the foamyHexMeshDict and"
        << " writing a re-sampled surface to " << exportName
        << nl << endl;

    cpuTime timer;

    IOdictionary foamyHexMeshDict
    (
        IOobject
        (
            "foamyHexMeshDict",
            runTime.system(),
            runTime,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    // Define/load all geometry
    searchableSurfaces allGeometry
    (
        IOobject
        (
            "cvSearchableSurfaces",
            runTime.constant(),
            "triSurface",
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        foamyHexMeshDict.subDict("geometry"),
        foamyHexMeshDict.lookupOrDefault("singleRegionName", true)
    );

    Info<< "Geometry read in = "
        << timer.cpuTimeIncrement() << " s." << nl << endl;


    Random rndGen(64293*Pstream::myProcNo());

    conformationSurfaces geometryToConformTo
    (
        runTime,
        rndGen,
        allGeometry,
        foamyHexMeshDict.subDict("surfaceConformation")
    );

    Info<< "Set up geometry in = "
        << timer.cpuTimeIncrement() << " s." << nl << endl;


    const searchableSurfaces& geometry = geometryToConformTo.geometry();
    const labelList& surfaces = geometryToConformTo.surfaces();


    const label minLevel = 2;

    // The max cube size follows from the minLevel and the default cube size
    // (1)
    const scalar maxSize = 1.0 / (1 << minLevel);
    const scalar halfMaxSize = 0.5*maxSize;


    // Scale the geometry to fit within
    // halfMaxSize .. 1-halfMaxSize

    scalar wantedRange = 1.0-maxSize;

    const treeBoundBox bb = geometryToConformTo.globalBounds();

    const vector scale = cmptDivide
    (
        point(wantedRange, wantedRange, wantedRange),
        bb.span()
    );
    const vector offset =
        cmptDivide
        (
            point(halfMaxSize, halfMaxSize, halfMaxSize),
            scale
        )
       -bb.min();


    const pointConversion converter(scale, offset);


    // Marching cubes

    OptOctree octree;

    distanceCalc ref
    (
        8,          // maxLevel
        0.0,        // distance
        minLevel,   // minLevel
        geometryToConformTo,
        converter
    );

    octree.refine(&ref);
    octree.set_impl(&ref);

    Info<< "Calculated octree in = "
        << timer.cpuTimeIncrement() << " s." << nl << endl;

    MarchingCubes& mc = octree.mc();

    mc.clean_all() ;
    octree.build_isosurface(&ref) ;

    Info<< "Constructed iso surface of distance in = "
        << timer.cpuTimeIncrement() << " s." << nl << endl;

    // Write output file
    if (mc.ntrigs() > 0)
    {
        Triangle* triangles = mc.triangles();
        label nTris = mc.ntrigs();
        Foam::DynamicList<labelledTri> tris(mc.ntrigs());
        for (label triI = 0; triI < nTris; ++triI)
        {
            const Triangle& t = triangles[triI];
            if (t.v1 != t.v2 && t.v1 != t.v3 && t.v2 != t.v3)
            {
                tris.append
                (
                    labelledTri
                    (
                        triangles[triI].v1,
                        triangles[triI].v2,
                        triangles[triI].v3,
                        0                       // region
                    )
                );
            }
        }


        Point* vertices = mc.vertices();
        pointField points(mc.nverts());
        forAll(points, pointi)
        {
            const Point& v = vertices[pointi];
            points[pointi] = converter.toGlobal(v);
        }


        // Find correspondence to original surfaces
        labelList regionOffsets(surfaces.size());
        label nRegions = 0;
        forAll(surfaces, i)
        {
            const wordList& regions = geometry[surfaces[i]].regions();
            regionOffsets[i] = nRegions;
            nRegions += regions.size();
        }


        geometricSurfacePatchList patches(nRegions);
        nRegions = 0;
        forAll(surfaces, i)
        {
            const wordList& regions = geometry[surfaces[i]].regions();

            forAll(regions, regionI)
            {
                patches[nRegions] = geometricSurfacePatch
                (
                    "patch",
                    geometry[surfaces[i]].name() + "_" + regions[regionI],
                    nRegions
                );
                nRegions++;
            }
        }

        triSurface s(tris, patches, points, true);
        tris.clearStorage();

        Info<< "Extracted triSurface in = "
            << timer.cpuTimeIncrement() << " s." << nl << endl;


        // Find out region on local surface of nearest point
        {
            Foam::List<pointIndexHit> hitInfo;
            labelList hitSurfaces;
            geometryToConformTo.findSurfaceNearest
            (
                s.faceCentres(),
                scalarField(s.size(), sqr(great)),
                hitInfo,
                hitSurfaces
            );

            // Get region
            DynamicList<pointIndexHit> surfInfo(hitSurfaces.size());
            DynamicList<label> surfIndices(hitSurfaces.size());

            forAll(surfaces, surfI)
            {
                // Extract info on this surface
                surfInfo.clear();
                surfIndices.clear();
                forAll(hitSurfaces, triI)
                {
                    if (hitSurfaces[triI] == surfI)
                    {
                        surfInfo.append(hitInfo[triI]);
                        surfIndices.append(triI);
                    }
                }

                // Calculate sideness of these surface points
                labelList region;
                geometry[surfaces[surfI]].getRegion(surfInfo, region);

                forAll(region, i)
                {
                    label triI = surfIndices[i];
                    s[triI].region() = regionOffsets[surfI]+region[i];
                }
            }
        }

        Info<< "Re-patched surface in = "
            << timer.cpuTimeIncrement() << " s." << nl << endl;

        triSurfaceMesh smesh
        (
            IOobject
            (
                exportName,
                runTime.constant(), // instance
                "triSurface",
                runTime,            // registry
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            s
        );

        Info<< "writing surfMesh:\n  "
            << smesh.searchableSurface::objectPath() << nl << endl;
        smesh.searchableSurface::write();

        Info<< "Written surface in = "
            << timer.cpuTimeIncrement() << " s." << nl << endl;
    }

    mc.clean_all() ;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
