singlePhaseTransportModel laminarTransport(U, phi);

// Thermal expansion coefficient [1/K]
dimensionedScalar beta
(
    "beta",
    dimless/dimTemperature,
    laminarTransport.lookup("beta")
);

// Reference temperature [K]
dimensionedScalar TRef("TRef", dimTemperature, laminarTransport.lookup("TRef"));

// Laminar Prandtl number
dimensionedScalar Pr("Pr", dimless, laminarTransport.lookup("Pr"));

// Turbulent Prandtl number
dimensionedScalar Prt("Prt", dimless, laminarTransport.lookup("Prt"));
