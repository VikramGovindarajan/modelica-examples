within Vessels;

partial model PartialLumpedVessel
"Lumped volume with a vector of fluid ports and replaceable heat transfer model"
  extends Interfaces.PartialLumpedVolume;
  
  import SI=Modelica.Units.SI;
  
  // Port definitions
  parameter Integer nPorts=0 "Number of ports";

  VesselFluidPorts_b ports[nPorts](redeclare each package Medium = Medium)
  "Fluid inlets and outlets";

  // Port properties
  parameter Boolean use_portsData=true
  "= false to neglect pressure loss and kinetic energy";
  parameter Vessels.BaseClasses.VesselPortsData[if use_portsData then nPorts else 0]
  portsData "Data of inlet/outlet ports";
  parameter Medium.MassFlowRate m_flow_nominal = if system.use_eps_Re then system.m_flow_nominal else 1e2*system.m_flow_small
  "Nominal value for mass flow rates in ports";
  parameter SI.MassFlowRate m_flow_small(min=0) = if system.use_eps_Re then system.eps_m_flow*m_flow_nominal else system.m_flow_small
  "Regularization range at zero mass flow rate";
  parameter Boolean use_Re = system.use_eps_Re
  "= true, if turbulent region is defined by Re, otherwise by m_flow_small";

  Medium.EnthalpyFlowRate ports_H_flow[nPorts];
  Medium.MassFlowRate ports_mXi_flow[nPorts,Medium.nXi];
  Medium.MassFlowRate[Medium.nXi] sum_ports_mXi_flow
  "Substance mass flows through ports";
  Medium.ExtraPropertyFlowRate ports_mC_flow[nPorts,Medium.nC];
  Medium.ExtraPropertyFlowRate[Medium.nC] sum_ports_mC_flow
  "Trace substance mass flows through ports";

  // Heat transfer through boundary
  parameter Boolean use_HeatTransfer = false
  "= true to use the HeatTransfer model";
  replaceable model HeatTransfer =
      Vessels.BaseClasses.HeatTransfer.IdealHeatTransfer
    constrainedby
  Vessels.BaseClasses.HeatTransfer.PartialVesselHeatTransfer
  "Wall heat transfer";
  HeatTransfer heatTransfer(
    redeclare package Medium = Medium,
    final n=1,
    final states = {medium.state},
    final use_k = use_HeatTransfer);
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort if use_HeatTransfer;

  // Conservation of kinetic energy
  Medium.Density[nPorts] portInDensities
  "Densities of the fluid at the device boundary";
  SI.Velocity[nPorts] portVelocities
  "Velocities of fluid flow at device boundary";
  SI.EnergyFlowRate[nPorts] ports_E_flow
  "Flow of kinetic and potential energy at device boundary";

  // Note: should use fluidLevel_start - portsData.height
  Real[nPorts] s(each start = fluidLevel_max)
  "Curve parameters for port flows vs. port pressures; for further details see, Modelica Tutorial: Ideal switching devices";
  Real[nPorts] ports_penetration
  "Penetration of port with fluid, depending on fluid level and port diameter";

  // treatment of pressure losses at ports
  SI.Area[nPorts] portAreas = {Modelica.Constants.pi/4*portsData_diameter[i]^2 for i in 1:nPorts};
  Medium.AbsolutePressure[nPorts] vessel_ps_static
  "Static pressures inside the vessel at the height of the corresponding ports, zero flow velocity";

  // determination of turbulent region
  constant SI.ReynoldsNumber Re_turbulent = 100 "cf. suddenExpansion";
  SI.MassFlowRate[nPorts] m_flow_turbulent;

protected
  input SI.Height fluidLevel = 0
  "Level of fluid in the vessel for treating heights of ports";
  parameter SI.Height fluidLevel_max = 1
  "Maximum level of fluid in the vessel";
  parameter SI.Area vesselArea = Modelica.Constants.inf
  "Area of the vessel used to relate to cross flow area of ports";

  // Treatment of use_portsData=false to neglect portsData and to not require its specification either in this case.
  // Remove portsData conditionally if use_portsData=false. Simplify their use in model equations by always
  // providing portsData_diameter and portsData_height, independent of the use_portsData setting.
  // Note: this moreover serves as work-around if a tool does not support a zero sized portsData record.
  Modelica.Blocks.Interfaces.RealInput[nPorts]
  portsData_diameter_internal = portsData.diameter if use_portsData and nPorts > 0;
  Modelica.Blocks.Interfaces.RealInput[nPorts] portsData_height_internal = portsData.height if use_portsData and nPorts > 0;
  Modelica.Blocks.Interfaces.RealInput[nPorts] portsData_zeta_in_internal = portsData.zeta_in if use_portsData and nPorts > 0;
  Modelica.Blocks.Interfaces.RealInput[nPorts] portsData_zeta_out_internal = portsData.zeta_out if use_portsData and nPorts > 0;
  Modelica.Blocks.Interfaces.RealInput[nPorts] portsData_diameter;
  Modelica.Blocks.Interfaces.RealInput[nPorts] portsData_height;
  Modelica.Blocks.Interfaces.RealInput[nPorts] portsData_zeta_in;
  Modelica.Blocks.Interfaces.RealInput[nPorts] portsData_zeta_out;
  Modelica.Blocks.Interfaces.BooleanInput[nPorts] regularFlow(each start=true);
  Modelica.Blocks.Interfaces.BooleanInput[nPorts] inFlow(each start=false);

equation
  mb_flow = sum(ports.m_flow);
  mbXi_flow = sum_ports_mXi_flow;
  mbC_flow  = sum_ports_mC_flow;
  Hb_flow = sum(ports_H_flow) + sum(ports_E_flow);
  Qb_flow = heatTransfer.Q_flows[1];

  // Only one connection allowed to a port to avoid unwanted ideal mixing
  for i in 1:nPorts loop
    assert(cardinality(ports[i]) <= 1,"
each ports[i] of volume can at most be connected to one component.
If two or more connections are present, ideal mixing takes
place with these connections, which is usually not the intention
of the modeller. Increase nPorts to add an additional port.
");
  end for;
  // Check for correct solution
  assert(fluidLevel <= fluidLevel_max, "Vessel is overflowing (fluidLevel > fluidLevel_max = " + String(fluidLevel) + ")");
  assert(fluidLevel > -1e-6*fluidLevel_max, "Fluid level (= " + String(fluidLevel) + ") is below zero meaning that the solution failed.");

  // Boundary conditions

  // treatment of conditional portsData
  connect(portsData_diameter, portsData_diameter_internal);
  connect(portsData_height, portsData_height_internal);
  connect(portsData_zeta_in, portsData_zeta_in_internal);
  connect(portsData_zeta_out, portsData_zeta_out_internal);
  if not use_portsData then
    portsData_diameter = zeros(nPorts);
    portsData_height = zeros(nPorts);
    portsData_zeta_in = zeros(nPorts);
    portsData_zeta_out = zeros(nPorts);
  end if;

  // actual definition of port variables
  for i in 1:nPorts loop
    portInDensities[i] = Medium.density(Medium.setState_phX(vessel_ps_static[i], inStream(ports[i].h_outflow), inStream(ports[i].Xi_outflow)));
    if use_portsData then
      // dp = 0.5*zeta*d*v*|v|
      // Note: assume vessel_ps_static for portVelocities to avoid algebraic loops for ports.p
      portVelocities[i] = smooth(0, ports[i].m_flow/portAreas[i]/Medium.density(Medium.setState_phX(vessel_ps_static[i], actualStream(ports[i].h_outflow), actualStream(ports[i].Xi_outflow))));
      // Note: the penetration should not go too close to zero as this would prevent a vessel from running empty
      ports_penetration[i] = Utilities.regStep(fluidLevel - portsData_height[i] - 0.1*portsData_diameter[i], 1, 1e-3, 0.1*portsData_diameter[i]);
      m_flow_turbulent[i]=if not use_Re then m_flow_small else
        max(m_flow_small, (Modelica.Constants.pi/8)*portsData_diameter[i]
                           *(Medium.dynamicViscosity(Medium.setState_phX(vessel_ps_static[i], inStream(ports[i].h_outflow), inStream(ports[i].Xi_outflow)))
                             + Medium.dynamicViscosity(medium.state))*Re_turbulent);
    else
      // an infinite port diameter is assumed
      portVelocities[i] = 0;
      ports_penetration[i] = 1;
      m_flow_turbulent[i] = Modelica.Constants.inf;
    end if;

    // fluid flow through ports
    regularFlow[i] = fluidLevel >= portsData_height[i];
    inFlow[i]      = not regularFlow[i] and (s[i] > 0 or portsData_height[i] >= fluidLevel_max);
    if regularFlow[i] then
      // regular operation: fluidLevel is above ports[i]
      // Note: >= covers default values of zero as well
      if use_portsData then
        /* Without regularization
           ports[i].p = vessel_ps_static[i] + 0.5*ports[i].m_flow^2/portAreas[i]^2
                        * noEvent(if ports[i].m_flow>0 then zeta_in[i]/portInDensities[i] else -zeta_out[i]/medium.d);
        */

        ports[i].p = vessel_ps_static[i] + (0.5/portAreas[i]^2*Utilities.regSquare2(ports[i].m_flow, m_flow_turbulent[i],
                          (portsData_zeta_in[i] - 1 + portAreas[i]^2/vesselArea^2)/portInDensities[i]*ports_penetration[i],
                          (portsData_zeta_out[i] + 1 - portAreas[i]^2/vesselArea^2)/medium.d/ports_penetration[i]));
        /*
          // alternative formulation m_flow=f(dp); not allowing the ideal portsData_zeta_in[i]=1 though
          ports[i].m_flow = smooth(2, portAreas[i]*Utilities.regRoot2(ports[i].p - vessel_ps_static[i], dp_small,
                                 2*portInDensities[i]/portsData_zeta_in[i],
                                 2*medium.d/portsData_zeta_out[i]));
        */
      else
        ports[i].p = vessel_ps_static[i];
      end if;
      s[i] = fluidLevel - portsData_height[i];

    elseif inFlow[i] then
      // ports[i] is above fluidLevel and has inflow
      ports[i].p = vessel_ps_static[i];
      s[i] = ports[i].m_flow;

    else
      // ports[i] is above fluidLevel, preventing outflow
      ports[i].m_flow = 0;
      s[i] = (ports[i].p - vessel_ps_static[i])/Medium.p_default*(portsData_height[i] - fluidLevel);
    end if;

    ports[i].h_outflow  = medium.h;
    ports[i].Xi_outflow = medium.Xi;
    ports[i].C_outflow  = C;

    ports_H_flow[i] = ports[i].m_flow * actualStream(ports[i].h_outflow)
    "Enthalpy flow";
    ports_E_flow[i] = ports[i].m_flow*(0.5*portVelocities[i]*portVelocities[i] + system.g*portsData_height[i])
    "Flow of kinetic and potential energy";
    ports_mXi_flow[i,:] = ports[i].m_flow * actualStream(ports[i].Xi_outflow)
    "Component mass flow";
    ports_mC_flow[i,:]  = ports[i].m_flow * actualStream(ports[i].C_outflow)
    "Trace substance mass flow";
  end for;

  for i in 1:Medium.nXi loop
    sum_ports_mXi_flow[i] = sum(ports_mXi_flow[:,i]);
  end for;

  for i in 1:Medium.nC loop
    sum_ports_mC_flow[i]  = sum(ports_mC_flow[:,i]);
  end for;

  connect(heatPort, heatTransfer.heatPorts[1]);
 
end PartialLumpedVessel;
