within Utilities;

function regStep
  "Approximation of a general step, such that the characteristic is continuous and differentiable"
  // extends Modelica.Icons.Function;
  input Real x "Abscissa value";
  input Real y1 "Ordinate value for x > 0";
  input Real y2 "Ordinate value for x < 0";
  input Real x_small(min=0) = 1e-5
    "Approximation of step for -x_small <= x <= x_small; x_small >= 0 required";
  output Real y "Ordinate value to approximate y = if x > 0 then y1 else y2";
algorithm
  y := smooth(1, if x >  x_small then y1 else
                 if x < -x_small then y2 else
                 if x_small > 0 then (x/x_small)*((x/x_small)^2 - 3)*(y2-y1)/4 + (y1+y2)/2 else (y1+y2)/2);
  
end regStep;
