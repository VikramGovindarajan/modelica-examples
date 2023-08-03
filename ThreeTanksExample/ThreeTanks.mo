model ThreeTanks "Demonstrating the usage of SimpleTank"

   replaceable package Medium =
      Modelica.Media.Water.ConstantPropertyLiquidWater                           constrainedby
    Modelica.Media.Interfaces.PartialMedium "Medium in the component";

  Vessels.OpenTank tank1(
    crossArea=1,
    redeclare package Medium = Medium,
    height=12,
    level_start=8,
    nPorts=1,
    portsData={Vessels.BaseClasses.VesselPortsData(diameter=
        0.1)});

  Vessels.OpenTank tank2(
    crossArea=1,
    redeclare package Medium = Medium,
    height=12,
    level_start=3,
    nPorts=1,
    portsData={Vessels.BaseClasses.VesselPortsData(diameter=
        0.1)});

  inner System system(energyDynamics=Types.Dynamics.FixedInitial);

  Vessels.OpenTank tank3(
    crossArea=1,
    redeclare package Medium = Medium,
    height=12,
    level_start=3,
    nPorts=1,
    portsData={Vessels.BaseClasses.VesselPortsData(diameter=
        0.1)});

  Pipes.StaticPipe pipe1(                    redeclare package
      Medium =                                                                       Medium,
    height_ab=2,
    length=2,
    diameter=0.1);

  Pipes.StaticPipe pipe2(                    redeclare package
      Medium =                                                                       Medium,
    height_ab=2,
    length=2,
    diameter=0.1);

  Pipes.StaticPipe pipe3(                    redeclare package
      Medium =                                                                       Medium,
    height_ab=-1,
    length=2,
    diameter=0.1);
equation
  connect(pipe1.port_a, pipe2.port_a);
  connect(pipe2.port_a, pipe3.port_a);
  connect(pipe3.port_b, tank3.ports[1]);
  connect(pipe1.port_b, tank1.ports[1]);
  connect(pipe2.port_b, tank2.ports[1]);

end ThreeTanks;
