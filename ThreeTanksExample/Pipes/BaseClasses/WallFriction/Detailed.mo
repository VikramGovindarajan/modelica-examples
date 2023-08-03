within Pipes.BaseClasses.WallFriction;

package Detailed
  "Pipe wall friction for laminar and turbulent flow (detailed characteristic)"

  extends PartialWallFriction(
            final use_mu = true,
            final use_roughness = true,
            final use_dp_small = true,
            final use_m_flow_small = true,
            final use_Re_turbulent = true);

  import ln = Modelica.Math.log "Logarithm, base e";
  import Modelica.Math.log10 "Logarithm, base 10";
  import Modelica.Math.exp "Exponential function";

  redeclare function extends massFlowRate_dp
    "Return mass flow rate m_flow as function of pressure loss dp, i.e., m_flow = f(dp), due to wall friction"
    import Modelica.Math;
	import SI=Modelica.Units.SI;
  protected
    Real Delta(min=0) = roughness/diameter "Relative roughness";
    SI.ReynoldsNumber Re1 = min((745*Math.exp(if Delta <= 0.0065 then 1 else 0.0065/Delta))^0.97, Re_turbulent)
      "Re leaving laminar curve";
    SI.ReynoldsNumber Re2 = Re_turbulent "Re entering turbulent curve";
    SI.DynamicViscosity mu "Upstream viscosity";
    SI.Density rho "Upstream density";
    SI.ReynoldsNumber Re "Reynolds number";
    Real lambda2 "Modified friction coefficient (= lambda*Re^2)";
    function interpolateInRegion2 = Modelica.Fluid.Dissipation.Utilities.Functions.General.CubicInterpolation_Re
      "Cubic Hermite spline interpolation in transition region";
  algorithm
    // Determine upstream density, upstream viscosity, and lambda2
    rho     := if dp >= 0 then rho_a else rho_b;
    mu      := if dp >= 0 then mu_a else mu_b;
    lambda2 := abs(dp)*2*diameter^3*rho/(length*mu*mu);

    // Determine Re under the assumption of laminar flow
    Re := lambda2/64;

    // Modify Re, if turbulent flow
    if Re > Re1 then
       Re :=-2*sqrt(lambda2)*Math.log10(2.51/sqrt(lambda2) + 0.27*Delta);
       if Re < Re2 then
          Re := interpolateInRegion2(Re, Re1, Re2, Delta, lambda2);
       end if;
    end if;

    // Determine mass flow rate
    m_flow := crossArea/diameter*mu*(if dp >= 0 then Re else -Re);
            
  end massFlowRate_dp;

  redeclare function extends pressureLoss_m_flow
    "Return pressure loss dp as function of mass flow rate m_flow, i.e., dp = f(m_flow), due to wall friction"
    import Modelica.Math;
    import Modelica.Constants.pi;
	import SI=Modelica.Units.SI;
  protected
    Real Delta(min=0) = roughness/diameter "Relative roughness";
    SI.ReynoldsNumber Re1 = min(745*Math.exp(if Delta <= 0.0065 then 1 else 0.0065/Delta), Re_turbulent)
      "Re leaving laminar curve";
    SI.ReynoldsNumber Re2 = Re_turbulent "Re entering turbulent curve";
    SI.DynamicViscosity mu "Upstream viscosity";
    SI.Density rho "Upstream density";
    SI.ReynoldsNumber Re "Reynolds number";
    Real lambda2 "Modified friction coefficient (= lambda*Re^2)";
    function interpolateInRegion2 = Modelica.Fluid.Dissipation.Utilities.Functions.General.CubicInterpolation_lambda
      "Cubic Hermite spline interpolation in transition region";
  algorithm
    // Determine upstream density and upstream viscosity
    rho     :=if m_flow >= 0 then rho_a else rho_b;
    mu      :=if m_flow >= 0 then mu_a else mu_b;

    // Determine Re, lambda2 and pressure drop
    Re := diameter*abs(m_flow)/(crossArea*mu);
    lambda2 := if Re <= Re1 then 64*Re else
              (if Re >= Re2 then 0.25*(Re/Math.log10(Delta/3.7 + 5.74/Re^0.9))^2 else
               interpolateInRegion2(Re, Re1, Re2, Delta));
    dp :=length*mu*mu/(2*rho*diameter*diameter*diameter)*
         (if m_flow >= 0 then lambda2 else -lambda2);
            
  end pressureLoss_m_flow;

  redeclare function extends massFlowRate_dp_staticHead
    "Return mass flow rate m_flow as function of pressure loss dp, i.e., m_flow = f(dp), due to wall friction and static head"
    import SI=Modelica.Units.SI;
	
  protected
    Real Delta(min=0) = roughness/diameter "Relative roughness";
    SI.ReynoldsNumber Re "Reynolds number";
    SI.ReynoldsNumber Re1 = min((745*exp(if Delta <= 0.0065 then 1 else 0.0065/Delta))^0.97, Re_turbulent)
      "Boundary between laminar regime and transition";
    SI.ReynoldsNumber Re2 = Re_turbulent
      "Boundary between transition and turbulent regime";
    SI.Pressure dp_a
      "Upper end of regularization domain of the m_flow(dp) relation";
    SI.Pressure dp_b
      "Lower end of regularization domain of the m_flow(dp) relation";
    SI.MassFlowRate m_flow_a
      "Value at upper end of regularization domain";
    SI.MassFlowRate m_flow_b
      "Value at lower end of regularization domain";

    SI.MassFlowRate dm_flow_ddp_fric_a
      "Derivative at upper end of regularization domain";
    SI.MassFlowRate dm_flow_ddp_fric_b
      "Derivative at lower end of regularization domain";

    SI.Pressure dp_grav_a = g_times_height_ab*rho_a
      "Static head if mass flows in design direction (a to b)";
    SI.Pressure dp_grav_b = g_times_height_ab*rho_b
      "Static head if mass flows against design direction (b to a)";

    // Properly define zero mass flow conditions
    SI.MassFlowRate m_flow_zero = 0;
    SI.Pressure dp_zero = (dp_grav_a + dp_grav_b)/2;
    Real dm_flow_ddp_fric_zero;

  algorithm
    dp_a := max(dp_grav_a, dp_grav_b)+dp_small;
    dp_b := min(dp_grav_a, dp_grav_b)-dp_small;

    if dp>=dp_a then
      // Positive flow outside regularization
      m_flow := Internal.m_flow_of_dp_fric(dp-dp_grav_a, rho_a, rho_b, mu_a, mu_b, length, diameter, crossArea, Re1, Re2, Delta);
    elseif dp<=dp_b then
      // Negative flow outside regularization
      m_flow := Internal.m_flow_of_dp_fric(dp-dp_grav_b, rho_a, rho_b, mu_a, mu_b, length, diameter, crossArea, Re1, Re2, Delta);
    else
      // Regularization parameters
      (m_flow_a, dm_flow_ddp_fric_a) := Internal.m_flow_of_dp_fric(dp_a-dp_grav_a, rho_a, rho_b, mu_a, mu_b, length, diameter, crossArea, Re1, Re2, Delta);
      (m_flow_b, dm_flow_ddp_fric_b) := Internal.m_flow_of_dp_fric(dp_b-dp_grav_b, rho_a, rho_b, mu_a, mu_b, length, diameter, crossArea, Re1, Re2, Delta);
      // Include a properly defined zero mass flow point
      // Obtain a suitable slope from the linear section slope c (value of m_flow is overwritten later)
      (m_flow, dm_flow_ddp_fric_zero) := Modelica.Fluid.Utilities.regFun3(dp_zero, dp_b, dp_a, m_flow_b, m_flow_a, dm_flow_ddp_fric_b, dm_flow_ddp_fric_a);
      // Do regularization
      if dp>dp_zero then
        m_flow := Modelica.Fluid.Utilities.regFun3(dp, dp_zero, dp_a, m_flow_zero, m_flow_a, dm_flow_ddp_fric_zero, dm_flow_ddp_fric_a);
      else
        m_flow := Modelica.Fluid.Utilities.regFun3(dp, dp_b, dp_zero, m_flow_b, m_flow_zero, dm_flow_ddp_fric_b, dm_flow_ddp_fric_zero);
      end if;
    end if;
    
  end massFlowRate_dp_staticHead;

  redeclare function extends pressureLoss_m_flow_staticHead
    "Return pressure loss dp as function of mass flow rate m_flow, i.e., dp = f(m_flow), due to wall friction and static head"
    import SI=Modelica.Units.SI;

  protected
    Real Delta(min=0) = roughness/diameter "Relative roughness";
    SI.ReynoldsNumber Re1 = min(745*exp(if Delta <= 0.0065 then 1 else 0.0065/Delta), Re_turbulent)
      "Boundary between laminar regime and transition";
    SI.ReynoldsNumber Re2 = Re_turbulent
      "Boundary between transition and turbulent regime";

    SI.MassFlowRate m_flow_a
      "Upper end of regularization domain of the dp(m_flow) relation";
    SI.MassFlowRate m_flow_b
      "Lower end of regularization domain of the dp(m_flow) relation";

    SI.Pressure dp_a "Value at upper end of regularization domain";
    SI.Pressure dp_b "Value at lower end of regularization domain";

    SI.Pressure dp_grav_a = g_times_height_ab*rho_a
      "Static head if mass flows in design direction (a to b)";
    SI.Pressure dp_grav_b = g_times_height_ab*rho_b
      "Static head if mass flows against design direction (b to a)";

    Real ddp_dm_flow_a
      "Derivative of pressure drop with mass flow rate at m_flow_a";
    Real ddp_dm_flow_b
      "Derivative of pressure drop with mass flow rate at m_flow_b";

    // Properly define zero mass flow conditions
    SI.MassFlowRate m_flow_zero = 0;
    SI.Pressure dp_zero = (dp_grav_a + dp_grav_b)/2;
    Real ddp_dm_flow_zero;

  algorithm
    m_flow_a := if dp_grav_a<dp_grav_b then
      Internal.m_flow_of_dp_fric(dp_grav_b - dp_grav_a, rho_a, rho_b, mu_a, mu_b, length, diameter, crossArea, Re1, Re2, Delta)+m_flow_small else
      m_flow_small;
    m_flow_b := if dp_grav_a<dp_grav_b then
      Internal.m_flow_of_dp_fric(dp_grav_a - dp_grav_b, rho_a, rho_b, mu_a, mu_b, length, diameter, crossArea, Re1, Re2, Delta)-m_flow_small else
      -m_flow_small;

    if m_flow>=m_flow_a then
      // Positive flow outside regularization
      dp := Internal.dp_fric_of_m_flow(m_flow, rho_a, rho_b, mu_a, mu_b, length, diameter, crossArea, Re1, Re2, Delta) + dp_grav_a;
    elseif m_flow<=m_flow_b then
      // Negative flow outside regularization
      dp := Internal.dp_fric_of_m_flow(m_flow, rho_a, rho_b, mu_a, mu_b, length, diameter, crossArea, Re1, Re2, Delta) + dp_grav_b;
    else
      // Regularization parameters
      (dp_a, ddp_dm_flow_a) := Internal.dp_fric_of_m_flow(m_flow_a, rho_a, rho_b, mu_a, mu_b, length, diameter, crossArea, Re1, Re2, Delta);
      dp_a := dp_a + dp_grav_a "Adding dp_grav to dp_fric to get dp";
      (dp_b, ddp_dm_flow_b) := Internal.dp_fric_of_m_flow(m_flow_b, rho_a, rho_b, mu_a, mu_b, length, diameter, crossArea, Re1, Re2, Delta);
      dp_b := dp_b + dp_grav_b "Adding dp_grav to dp_fric to get dp";
      // Include a properly defined zero mass flow point
      // Obtain a suitable slope from the linear section slope c (value of dp is overwritten later)
      (dp, ddp_dm_flow_zero) := Modelica.Fluid.Utilities.regFun3(m_flow_zero, m_flow_b, m_flow_a, dp_b, dp_a, ddp_dm_flow_b, ddp_dm_flow_a);
      // Do regularization
      if m_flow>m_flow_zero then
        dp := Modelica.Fluid.Utilities.regFun3(m_flow, m_flow_zero, m_flow_a, dp_zero, dp_a, ddp_dm_flow_zero, ddp_dm_flow_a);
      else
        dp := Modelica.Fluid.Utilities.regFun3(m_flow, m_flow_b, m_flow_zero, dp_b, dp_zero, ddp_dm_flow_b, ddp_dm_flow_zero);
      end if;
    end if;
    
  end pressureLoss_m_flow_staticHead;

package Internal
    "Functions to calculate mass flow rate from friction pressure drop and vice versa"
  extends Modelica.Icons.InternalPackage;
  function m_flow_of_dp_fric
      "Calculate mass flow rate as function of pressure drop due to friction"
    extends Modelica.Icons.Function;
    import SI=Modelica.Units.SI;
	
    input SI.Pressure dp_fric
        "Pressure loss due to friction (dp = port_a.p - port_b.p)";
    input SI.Density rho_a "Density at port_a";
    input SI.Density rho_b "Density at port_b";
    input SI.DynamicViscosity mu_a
        "Dynamic viscosity at port_a (dummy if use_mu = false)";
    input SI.DynamicViscosity mu_b
        "Dynamic viscosity at port_b (dummy if use_mu = false)";
    input SI.Length length "Length of pipe";
    input SI.Diameter diameter "Inner (hydraulic) diameter of pipe";
    input SI.Area crossArea "Inner cross section area";
    input SI.ReynoldsNumber Re1
        "Boundary between laminar regime and transition";
    input SI.ReynoldsNumber Re2
        "Boundary between transition and turbulent regime";
    input Real Delta(min=0) "Relative roughness";
    output SI.MassFlowRate m_flow "Mass flow rate from port_a to port_b";
    output Real dm_flow_ddp_fric
        "Derivative of mass flow rate with dp_fric";

    protected
    function interpolateInRegion2_withDerivative
        "Interpolation in log-log space using a cubic Hermite polynomial, where x=log10(lambda2), y=log10(Re)"
      extends Modelica.Icons.Function;

      input Real lambda2 "Known independent variable";
      input SI.ReynoldsNumber Re1
          "Boundary between laminar regime and transition";
      input SI.ReynoldsNumber Re2
          "Boundary between transition and turbulent regime";
      input Real Delta(min=0) "Relative roughness";
      input SI.Pressure dp_fric
          "Pressure loss due to friction (dp = port_a.p - port_b.p)";
      output SI.ReynoldsNumber Re "Unknown return variable";
      output Real dRe_ddp "Derivative of return value";
      // point lg(lambda2(Re1)) with derivative at lg(Re1)
      protected
      Real x1=log10(64*Re1);
      Real y1=log10(Re1);
      Real y1d=1;

      // Point lg(lambda2(Re2)) with derivative at lg(Re2)
      Real aux2=Delta/3.7 + 5.74/Re2^0.9;
      Real aux3=log10(aux2);
      Real L2=0.25*(Re2/aux3)^2;
      Real aux4=2.51/sqrt(L2) + 0.27*Delta;
      Real aux5=-2*sqrt(L2)*log10(aux4);
      Real x2=log10(L2);
      Real y2=log10(aux5);
      Real y2d=0.5 + (2.51/log(10))/(aux5*aux4);

      // Point of interest in transformed space
      Real x=log10(lambda2);
      Real y;
      Real dy_dx "Derivative in transformed space";
    algorithm
      // Interpolation
      (y, dy_dx) := Modelica.Fluid.Utilities.cubicHermite_withDerivative(x, x1, x2, y1, y2, y1d, y2d);

      // Return value
      Re := 10^y;

      // Derivative of return value
      dRe_ddp := Re/abs(dp_fric)*dy_dx;
      
    end interpolateInRegion2_withDerivative;

    SI.DynamicViscosity mu "Upstream viscosity";
    SI.Density rho "Upstream density";
    Real lambda2 "Modified friction coefficient (= lambda*Re^2)";
    SI.ReynoldsNumber Re "Reynolds number";
    Real dRe_ddp "dRe/ddp";
    Real aux1;
    Real aux2;

  algorithm
    // Determine upstream density and upstream viscosity
    if dp_fric >= 0 then
      rho := rho_a;
      mu  := mu_a;
    else
      rho := rho_b;
      mu  := mu_b;
    end if;

    // Positive mass flow rate
    lambda2 := abs(dp_fric)*2*diameter^3*rho/(length*mu*mu)
        "Known as lambda2=f(dp)";

    aux1:=(2*diameter^3*rho)/(length*mu^2);

    // Determine Re and dRe/ddp under the assumption of laminar flow
    Re := lambda2/64 "Hagen-Poiseuille";
    dRe_ddp := aux1/64 "Hagen-Poiseuille";

    // Modify Re, if turbulent flow
    if Re > Re1 then
      Re :=-2*sqrt(lambda2)*log10(2.51/sqrt(lambda2) + 0.27*Delta)
          "Colebrook-White";
      aux2 := sqrt(aux1*abs(dp_fric));
      dRe_ddp := 1/log(10)*(-2*log(2.51/aux2+0.27*Delta)*aux1/(2*aux2)+2*2.51/(2*abs(dp_fric)*(2.51/aux2+0.27*Delta)));
      if Re < Re2 then
        (Re, dRe_ddp) := interpolateInRegion2_withDerivative(lambda2, Re1, Re2, Delta, dp_fric);
      end if;
    end if;

    // Determine mass flow rate
    m_flow := crossArea/diameter*mu*(if dp_fric >= 0 then Re else -Re);
    // Determine derivative of mass flow rate with dp_fric
    dm_flow_ddp_fric := crossArea/diameter*mu*dRe_ddp;
    
  end m_flow_of_dp_fric;

  function dp_fric_of_m_flow
      "Calculate pressure drop due to friction as function of mass flow rate"
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
    input SI.Area crossArea "Inner cross section area";
    input SI.ReynoldsNumber Re1
        "Boundary between laminar regime and transition";
    input SI.ReynoldsNumber Re2
        "Boundary between transition and turbulent regime";
    input Real Delta(min=0) "Relative roughness";
    output SI.Pressure dp_fric
        "Pressure loss due to friction (dp_fric = port_a.p - port_b.p - dp_grav)";
    output Real ddp_fric_dm_flow
        "Derivative of pressure drop with mass flow rate";

    protected
    function interpolateInRegion2
        "Interpolation in log-log space using a cubic Hermite polynomial, where x=log10(Re), y=log10(lambda2)"
      extends Modelica.Icons.Function;

      input SI.ReynoldsNumber Re "Known independent variable";
      input SI.ReynoldsNumber Re1
          "Boundary between laminar regime and transition";
      input SI.ReynoldsNumber Re2
          "Boundary between transition and turbulent regime";
      input Real Delta(min=0) "Relative roughness";
      input SI.MassFlowRate m_flow "Mass flow rate from port_a to port_b";
      output Real lambda2 "Unknown return value";
      output Real dlambda2_dm_flow "Derivative of return value";
      // point lg(lambda2(Re1)) with derivative at lg(Re1)
      protected
      Real x1 = log10(Re1);
      Real y1 = log10(64*Re1);
      Real y1d = 1;

      // Point lg(lambda2(Re2)) with derivative at lg(Re2)
      Real aux2 = Delta/3.7 + 5.74/Re2^0.9;
      Real aux3 = log10(aux2);
      Real L2 = 0.25*(Re2/aux3)^2;
      Real x2 = log10(Re2);
      Real y2 = log10(L2);
      Real y2d = 2+(2*5.74*0.9)/(log(aux2)*Re2^0.9*aux2);

      // Point of interest in transformed space
      Real x=log10(Re);
      Real y;
      Real dy_dx "Derivative in transformed space";
    algorithm
      // Interpolation
      (y, dy_dx) := Modelica.Fluid.Utilities.cubicHermite_withDerivative(x, x1, x2, y1, y2, y1d, y2d);

      // Return value
      lambda2 := 10^y;

      // Derivative of return value
      dlambda2_dm_flow := lambda2/abs(m_flow)*dy_dx;
      
    end interpolateInRegion2;

    SI.DynamicViscosity mu "Upstream viscosity";
    SI.Density rho "Upstream density";
    SI.ReynoldsNumber Re "Reynolds number";
    Real lambda2 "Modified friction coefficient (= lambda*Re^2)";
    Real dlambda2_dm_flow "dlambda2/dm_flow";
    Real aux1;
    Real aux2;

  algorithm
    // Determine upstream density and upstream viscosity
    if m_flow >= 0 then
      rho := rho_a;
      mu  := mu_a;
    else
      rho := rho_b;
      mu  := mu_b;
    end if;

    // Determine Reynolds number
    Re := abs(m_flow)*diameter/(crossArea*mu);

    aux1 := diameter/(crossArea*mu);

    // Use correlation for lambda2 depending on actual conditions
    if Re <= Re1 then
      lambda2 := 64*Re "Hagen-Poiseuille";
      dlambda2_dm_flow := 64*aux1 "Hagen-Poiseuille";
    elseif Re >= Re2 then
      lambda2 := 0.25*(Re/log10(Delta/3.7 + 5.74/Re^0.9))^2 "Swamee-Jain";
      aux2 := Delta/3.7+5.74/((aux1*abs(m_flow))^0.9);
      dlambda2_dm_flow := 0.5*aux1*Re*log(10)^2*(1/(log(aux2)^2)+(5.74*0.9)/(log(aux2)^3*Re^0.9*aux2))
          "Swamee-Jain";
    else
      (lambda2, dlambda2_dm_flow) := interpolateInRegion2(Re, Re1, Re2, Delta, m_flow);
    end if;

    // Compute pressure drop from lambda2
    dp_fric :=length*mu*mu/(2*rho*diameter*diameter*diameter)*
         (if m_flow >= 0 then lambda2 else -lambda2);

    // Compute derivative from dlambda2/dm_flow
    ddp_fric_dm_flow := (length*mu^2)/(2*diameter^3*rho)*dlambda2_dm_flow;
    
  end dp_fric_of_m_flow;
end Internal;

end Detailed;
