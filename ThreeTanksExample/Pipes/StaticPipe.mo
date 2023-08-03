within Pipes;

model StaticPipe "Basic pipe flow model without storage of mass or energy"

  // extending PartialStraightPipe
  extends Pipes.BaseClasses.PartialStraightPipe;

  // Initialization
  parameter Medium.AbsolutePressure p_a_start=system.p_start
    "Start value of pressure at port a";
  parameter Medium.AbsolutePressure p_b_start=p_a_start
    "Start value of pressure at port b";
  parameter Medium.MassFlowRate m_flow_start = system.m_flow_start
    "Start value for mass flow rate";

  FlowModel flowModel(
          redeclare package Medium = Medium,
          final n=2,
          states={Medium.setState_phX(port_a.p, inStream(port_a.h_outflow), inStream(port_a.Xi_outflow)),
                 Medium.setState_phX(port_b.p, inStream(port_b.h_outflow), inStream(port_b.Xi_outflow))},
          vs={port_a.m_flow/Medium.density(flowModel.states[1])/flowModel.crossAreas[1],
              -port_b.m_flow/Medium.density(flowModel.states[2])/flowModel.crossAreas[2]}/nParallel,
          final momentumDynamics=Types.Dynamics.SteadyState,
          final allowFlowReversal=allowFlowReversal,
          final p_a_start=p_a_start,
          final p_b_start=p_b_start,
          final m_flow_start=m_flow_start,
          final nParallel=nParallel,
          final pathLengths={length},
          final crossAreas={crossArea, crossArea},
          final dimensions={4*crossArea/perimeter, 4*crossArea/perimeter},
          final roughnesses={roughness, roughness},
          final dheights={height_ab},
          final g=system.g) "Flow model";
equation
  // Mass balance
  port_a.m_flow = flowModel.m_flows[1];
  0 = port_a.m_flow + port_b.m_flow;
  port_a.Xi_outflow = inStream(port_b.Xi_outflow);
  port_b.Xi_outflow = inStream(port_a.Xi_outflow);
  port_a.C_outflow = inStream(port_b.C_outflow);
  port_b.C_outflow = inStream(port_a.C_outflow);

  // Energy balance, considering change of potential energy
  // Wb_flow = v*A*dpdx + v*F_fric
  //         = m_flow/d/A * (A*dpdx + A*pressureLoss.dp_fg - F_grav)
  //         = m_flow/d/A * (-A*g*height_ab*d)
  //         = -m_flow*g*height_ab
  port_b.h_outflow = inStream(port_a.h_outflow) - system.g*height_ab;
  port_a.h_outflow = inStream(port_b.h_outflow) + system.g*height_ab;

end StaticPipe;
