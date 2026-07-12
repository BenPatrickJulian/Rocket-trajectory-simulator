# Rocket-trajectory-simulator
A custom rocket flight dynamics simulator written in python.

To use this code the user must have numpy, pandas and matplotlib.pyplot downloaded as well as a csv file containing the raw data of the thrust curve of the engine in the same directory and a csv file of openrockets raw data of the flight in the same directory if you want the comparison between the 2. The variables at the top of the file indicated to be changed must be changed to the dimensions and qualities of the rocket being modelled.

## Features 
* Euler-cromer integration is used as a first order numerical method to determine position from acceleration.;
* customisable rocket: Many different parts of the rocket are changeable due to what rocket is being used.
* 1D plotting capabilities: Choice of altitude, velocity and acceleration can be plotted against time and the max of each of these shown. As well as this the curves can be plotted against openrocket and a normalized root mean square error is used to analyze the accuracy consistently being >95% in each of these curves.
* A variable air density calculator using the barometric formula.

## Physics engine
The simulator updates the acceleration every time period dt by taking into account gravity, drag and thrust from the engine.

$$\Sigma F = F_{\text{thrust}} - F_{\text{drag}} - F_{\text{gravity}}$$

Drag is calculated via the quadratic drag formula and the air density calculated by the barometric formula. Thrust is found by interpolating the thrust curve found at https://www.thrustcurve.org/ for the given engine. This acceleration is then numerically integrated to a velocity and position which is then recorded and plotted against time at the end. 

After consideration and research doing a component breakdown approach for drag coefficient will not be accurate enough to implement until the model uses 6 degrees of freedom therefore the improvement is moved until after updated to 6 degrees of freedom.

## Improvements to come 
- [ ] Runge-kutta 4th order integration is to be implemented over euler-cromer integration making it far more accurate and precise.
- [ ] The model is to be updated to 3 degrees of freedom to be able to model more advanced situations
- [ ] The model is to be updated to 6 degrees of freedom as to be more accurate and able to analyze stability as well as trajectory
- [ ] Drag coefficient is to be calculated using a component breakdown approach with 3 main points the skin friction drag, base drag and pressure drag allowing it to change dynamically throughout.
