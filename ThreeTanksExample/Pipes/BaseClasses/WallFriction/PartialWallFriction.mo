within Pipes.BaseClasses.WallFriction;

partial package PartialWallFriction
  "Partial wall friction characteristic (base package of all wall friction characteristics)"
  extends Modelica.Icons.Package;
  import Modelica.Constants.pi;

// Constants to be set in subpackages
  constant Boolean use_mu = true
    "= true, if mu_a/mu_b are used in function, otherwise value is not used";
  constant Boolean use_roughness = true
    "= true, if roughness is used in function, otherwise value is not used";
  constant Boolean use_dp_small = true
    "= true, if dp_small is used in function, otherwise value is not used";
  constant Boolean use_m_flow_small = true
    "= true, if m_flow_small is used in function, otherwise value is not used";
  constant Boolean dp_is_zero = false
    "= true, if no wall friction is present, i.e., dp = 0 (function massFlowRate_dp() cannot be used)";
  constant Boolean use_Re_turbulent = true
    "= true, if Re_turbulent input is used in function, otherwise value is not used";

// pressure loss characteristic functions
  replaceable partial function massFlowRate_dp
    "Return mass flow rate m_flow as function of pressure loss dp, i.e., m_flow = f(dp), due to wall friction"
    extends Modelica.Icons.Function;
	import SI=Modelica.Units.SI;

    input SI.Pressure dp "Pressure loss (dp = port_a.p - port_b.p)";
    input SI.Density rho_a "Density at port_a";
    input SI.Density rho_b "Density at port_b";
    input SI.DynamicViscosity mu_a
      "Dynamic viscosity at port_a (dummy if use_mu = false)";
    input SI.DynamicViscosity mu_b
      "Dynamic viscosity at port_b (dummy if use_mu = false)";
    input SI.Length length "Length of pipe";
    input SI.Diameter diameter "Inner (hydraulic) diameter of pipe";
    input SI.Area crossArea = pi*diameter^2/4 "Inner cross section area";
    input Types.Roughness roughness = 2.5e-5
      "Absolute roughness of pipe, with a default for a smooth steel pipe (dummy if use_roughness = false)";
    input SI.AbsolutePressure dp_small = 1
      "Regularization of zero flow if |dp| < dp_small (dummy if use_dp_small = false)";
    input SI.ReynoldsNumber Re_turbulent = 4000
      "Turbulent flow if Re >= Re_turbulent (dummy if use_Re_turbulent = false)";

    output SI.MassFlowRate m_flow "Mass flow rate from port_a to port_b";
  
  end massFlowRate_dp;

  replaceable partial function massFlowRate_dp_staticHead
    "Return mass flow rate m_flow as function of pressure loss dp, i.e., m_flow = f(dp), due to wall friction and static head"
    extends Modelica.Icons.Function;
    import SI=Modelica.Units.SI;
	
    input SI.Pressure dp "Pressure loss (dp = port_a.p - port_b.p)";
    input SI.Density rho_a "Density at port_a";
    input SI.Density rho_b "Density at port_b";
    input SI.DynamicViscosity mu_a
      "Dynamic viscosity at port_a (dummy if use_mu = false)";
    input SI.DynamicViscosity mu_b
      "Dynamic viscosity at port_b (dummy if use_mu = false)";
    input SI.Length length "Length of pipe";
    input SI.Diameter diameter "Inner (hydraulic) diameter of pipe";
    input Real g_times_height_ab(unit="m2/s2")
      "Gravity times (Height(port_b) - Height(port_a))";
    input SI.Area crossArea = pi*diameter^2/4 "Inner cross section area";
    input Types.Roughness roughness = 2.5e-5
      "Absolute roughness of pipe, with a default for a smooth steel pipe (dummy if use_roughness = false)";
    input SI.AbsolutePressure dp_small=1
      "Regularization of zero flow if |dp| < dp_small (dummy if use_dp_small = false)";
    input SI.ReynoldsNumber Re_turbulent = 4000
      "Turbulent flow if Re >= Re_turbulent (dummy if use_Re_turbulent = false)";

    output SI.MassFlowRate m_flow "Mass flow rate from port_a to port_b";
    
  end massFlowRate_dp_staticHead;

  replaceable partial function pressureLoss_m_flow
    "Return pressure loss dp as function of mass flow rate m_flow, i.e., dp = f(m_flow), due to wall friction"
    extends Modelica.Icons.Function;
    import SI=Modelica.Units.SI;
	
    input SI.MassFlowRate m_flow "Mass flow rate from port_a to port_b";
    input SI.Density rho_a "Density at port_a";
    input SI.Density rho_b "Density at port_b";
    input SI.DynamicViscosity mu_a
      "Dynamic viscosity at port_a (dummy if use_mu = false)";
    input SI.DynamicViscosity mu_b
      "Dynamic viscosity at port_b (dummy if use_mu = false)";
    input SI.Length length "Length of pipe";
    input SI.Diameter diameter "Inner (hydraulic) diameter of pipe";
    input SI.Area crossArea = pi*diameter^2/4 "Inner cross section area";
    input Types.Roughness roughness = 2.5e-5
      "Absolute roughness of pipe, with a default for a smooth steel pipe (dummy if use_roughness = false)";
    input SI.MassFlowRate m_flow_small = 0.01
      "Regularization of zero flow if |m_flow| < m_flow_small (dummy if use_m_flow_small = false)";
    input SI.ReynoldsNumber Re_turbulent = 4000
      "Turbulent flow if Re >= Re_turbulent (dummy if use_Re_turbulent = false)";

    output SI.Pressure dp "Pressure loss (dp = port_a.p - port_b.p)";

  
  end pressureLoss_m_flow;

  replaceable partial function pressureLoss_m_flow_staticHead
    "Return pressure loss dp as function of mass flow rate m_flow, i.e., dp = f(m_flow), due to wall friction and static head"
            extends Modelica.Icons.Function;
    import SI=Modelica.Units.SI;
	
	input SI.MassFlowRate m_flow "Mass flow rate from port_a to port_b";
    input SI.Density rho_a "Density at port_a";
    input SI.Density rho_b "Density at port_b";
    input SI.DynamicViscosity mu_a
      "Dynamic viscosity at port_a (dummy if use_mu = false)";
    input SI.DynamicViscosity mu_b
      "Dynamic viscosity at port_b (dummy if use_mu = false)";
    input SI.Length length "Length of pipe";
    input SI.Diameter diameter "Inner (hydraulic) diameter of pipe";
    input Real g_times_height_ab(unit="m2/s2")
      "Gravity times (Height(port_b) - Height(port_a))";
    input SI.Area crossArea = pi*diameter^2/4 "Inner cross section area";
    input Types.Roughness roughness = 2.5e-5
      "Absolute roughness of pipe, with a default for a smooth steel pipe (dummy if use_roughness = false)";
    input SI.MassFlowRate m_flow_small = 0.01
      "Regularization of zero flow if |m_flow| < m_flow_small (dummy if use_m_flow_small = false)";
    input SI.ReynoldsNumber Re_turbulent = 4000
      "Turbulent flow if Re >= Re_turbulent (dummy if use_Re_turbulent = false)";

    output SI.Pressure dp "Pressure loss (dp = port_a.p - port_b.p)";

  
  end pressureLoss_m_flow_staticHead;
  
end PartialWallFriction;
