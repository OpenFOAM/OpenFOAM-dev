Info<< "\nConstructing coal cloud" << endl;
coalCloud coalParcels
(
    "coalCloud1",
    rho,
    U,
    g,
    slgThermo
);

Info<< "\nConstructing limestone cloud" << endl;
basicThermoCloud limestoneParcels
(
    "limestoneCloud1",
    rho,
    U,
    g,
    slgThermo
);
