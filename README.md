# OPNA_2 - Sturm's Theorem Evaluator
```
./opna_2 3x^2+3x-4 -5 5 --verbose

 Sturm sequence    
-------------------
 P0(x) = 3x^2+3x-4 
 P1(x) = 6x+3      
 P2(x) = 19/4      

 Point   Sign sequence   Sign variations 
------- --------------- -----------------
 -5      (+,-,+)         2               
  5      (+,+,+)         0               

 Number of real roots of p0(x) on interval (-5,5] 
--------------------------------------------------
 V(-5) - V(5) = 2                                 
```
```
Usage: ./opna_2 [OPTIONS] <POLYNOMIAL> <INTERVAL_FROM> <INTERVAL_TO> 
Allowed options:
  --help                  print help message
  --version               print version information
  --examples              show examples
  --gcd                   apply on P(x)/Q(x) where Q(x)=GCD(P(x),P'(x))
  --verbose               detailed output of every step
  --table_style arg       Specify table style for verbose output: 
                          nice|double|simple|empty
  --polynomial_string arg polynomial_string to be parsed
  --from arg              interval from
  --to arg                interval to
```
### Examples
```
    Command                              Description                                                      
 1  ./opna_2 x^4+x^3-x-1 -5 5            Evaluate P(x) = x^4+x^3-x-1 with Sturm's theorem on range (-5,5] 
 2  ./opna_2 x^4+x^3-x-1 -5 5 --verbose  Same as above but with detailed output of every step             
 3  ./opna_2 1/2x^2-5.0                  Fractions and floats are supported as well                       
 4  ./opna_2 x^3+x^5+2x^3                Mixed ordering and repeating of the same element is allowed   
```
### Build requirements
 C++17 compatible compiler, boost libraries, cmake
