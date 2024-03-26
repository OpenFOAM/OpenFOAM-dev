Notes for fluentMeshToFoam with zone preservation
#################################################

1. New option added:
    - writeSets:
    Writes all Fluent boundaries faceSets preserving Fluent names
    Writes all Fluent regions to cellSets preserving Fluent names
    lines: 1375 - 1393 & 1673 - 1741
    sets are useful for post-processing using foamToVTK with the "-faceSet
    <name>" and "-cellSet <name>" options.

    - writeZones:
    Writes all regions to cellZones preserving Fluent names
    Writes all region internal face to faceZones preserving Fluent names
    lines: 1545 - 1667
    Zones are useful for porous media and MRF calculations

2. Zone Access
    - Zones are simple lists of label lists that can be accessed from polyMesh
    with the cellZones(), faceZones() and pointZones() member functions

    - Example (Members from polyMesh.H and Zones.H):
    const labelList& thisCellZone = mesh.cellZones()["thisZoneName"];

    - Zone integrity is preserved during mesh modification and decompomposition.

    - Once created via addZones, zones allow modification through non-const
    access

3. Fluent boundary types.
    - All internal and baffle elements are ignored during conversion

    - Boundary faces labelled as internal (i.e. interior, interface, internal,
    solid, fan, radiator, porous-jump) but that are in fact external boundaries
    will be added to a default wall boundary.
