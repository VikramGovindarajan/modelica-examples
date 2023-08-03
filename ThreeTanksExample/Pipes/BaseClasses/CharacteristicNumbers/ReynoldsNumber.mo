within Pipes.BaseClasses.CharacteristicNumbers;

function ReynoldsNumber "Return Reynolds number from v, rho, mu, D"
  // extends Modelica.Icons.Function;
  import SI=Modelica.Units.SI;
  
  input SI.Velocity v "Mean velocity of fluid flow";
  input SI.Density rho "Fluid density";
  input SI.DynamicViscosity mu "Dynamic (absolute) viscosity";
  input SI.Length D
    "Characteristic dimension (hydraulic diameter of pipes)";
  output SI.ReynoldsNumber Re "Reynolds number";
algorithm
  Re := abs(v)*rho*D/mu;

end ReynoldsNumber;
