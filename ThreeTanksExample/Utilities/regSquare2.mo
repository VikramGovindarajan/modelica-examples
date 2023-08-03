within Utilities;

function regSquare2
  "Anti-symmetric approximation of square with discontinuous factor so that the first derivative is non-zero and is continuous"
  // extends Modelica.Icons.Function;
  input Real x "Abscissa value";
  input Real x_small(min=0)=0.01
    "Approximation of function for |x| <= x_small";
  input Real k1(min=0)=1 "y = (if x>=0 then k1 else k2)*x*|x|";
  input Real k2(min=0)=1 "y = (if x>=0 then k1 else k2)*x*|x|";
  input Boolean use_yd0 = false "= true, if yd0 shall be used";
  input Real yd0(min=0)=1 "Desired derivative at x=0: dy/dx = yd0";
  output Real y "Ordinate value";
protected
  encapsulated function regSquare2_utility
    "Interpolating with two 3-order polynomials with a prescribed derivative at x=0"
    // import Modelica;
    // extends Modelica.Icons.Function;
    import Modelica.Fluid.Utilities.evaluatePoly3_derivativeAtZero;
     input Real x;
     input Real x1 "Approximation of function abs(x) < x1";
     input Real k1 "y = (if x>=0 then k1 else -k2)*x*|x|; k1 >= k2";
     input Real k2 "y = (if x>=0 then k1 else -k2)*x*|x|";
     input Boolean use_yd0 = false "= true, if yd0 shall be used";
     input Real yd0(min=0)=1 "Desired derivative at x=0: dy/dx = yd0";
     output Real y;
  protected
     Real x2;
     Real y1;
     Real y2;
     Real y1d;
     Real y2d;
     Real w;
     Real w1;
     Real w2;
     Real y0d;
     Real ww;
  algorithm
     // x2 :=-x1*(k2/k1)^2;
     x2 := -x1;
     if x <= x2 then
        y := -k2*x^2;
     else
         y1 := k1*x1^2;
         y2 :=-k2*x2^2;
        y1d := k1*2*x1;
        y2d :=-k2*2*x2;
        if use_yd0 then
           y0d :=yd0;
        else
           /* Determine derivative, such that first and second derivative
            of left and right polynomial are identical at x=0:
            see derivation in function regRoot2
         */
           w :=x2/x1;
           y0d := ( (3*y2 - x2*y2d)/w - (3*y1 - x1*y1d)*w) /(2*x1*(1 - w));
        end if;

        /* Modify derivative y0d, such that the polynomial is
         monotonically increasing. A sufficient condition is
           0 <= y0d <= sqrt(5)*k_i*|x_i|
      */
        w1 :=sqrt(5)*k1*x1;
        w2 :=sqrt(5)*k2*abs(x2);
        // y0d :=min(y0d, 0.9*min(w1, w2));
        ww :=0.9*(if w1 < w2 then w1 else w2);
        if ww < y0d then
           y0d :=ww;
        end if;
        y := if x >= 0 then evaluatePoly3_derivativeAtZero(x,x1,y1,y1d,y0d) else
                            evaluatePoly3_derivativeAtZero(x,x2,y2,y2d,y0d);
     end if;
  end regSquare2_utility;
algorithm
  y := smooth(2,if x >= x_small then k1*x^2 else
                if x <= -x_small then -k2*x^2 else
                if k1 >= k2 then regSquare2_utility(x,x_small,k1,k2,use_yd0,yd0) else
                                -regSquare2_utility(-x,x_small,k2,k1,use_yd0,yd0));
  
end regSquare2;
