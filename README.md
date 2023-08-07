# tennis-service

Tennis Service Trajectory Simulation
This MATLAB project aims to simulate the trajectory of a tennis service. The simulation takes the initial speed and direction of the tennis ball as input data and determines its position as a function of time during the service motion. The simulation adheres to classical tennis rules to determine whether the service is in or out of the court.

Tennis Court and Dimensions
The court boundaries and net function are used to ensure that the ball remains within the playable part of the court during the simulation.

Numerical Method: Adams-Moulton Method of Order 4
To model the trajectory of the ball, we utilize the Adams-Moulton method of order 4 to solve the differential equations governing the motion of the ball. This numerical method is compared with the analytical solution to validate the simulation's accuracy.

Inclusion of Physical Forces
Two scenarios are considered for the ball's motion:

Gravity-Only Motion: In this case, only the force of gravity is considered to determine the ball's trajectory.
Full Model: The simulation takes into account the combined effects of gravity, Magnus effect, and drift on the ball's trajectory.

Key Observations
Based on different data sets, the following observations were made:

The graphs of the analytical solution and the numerical solution coincide, validating the accuracy of the simulation.
The ball can clear the net either with a low initial speed and an angle of less than 90° or with a high initial speed and an angle of more than 90°.
The inclusion of the Magnus effect and drift in the simulation results in a slower ball movement.
How to Use
Set the initial speed and direction of the ball in the MATLAB script.
Run the script to perform the simulation and determine the trajectory of the ball.
The simulation considers the ball's motion until it goes out of bounds or hits the ground after the rebound.

Notes
This project is designed for educational and illustrative purposes, using simplified assumptions and models.
The Adams-Moulton method is used as a numerical solution, and the accuracy of the simulation depends on the input parameters and model limitations.

Author
The MATLAB code and README are authored by Polina Shamokhina, p.shamokhina@mail.ru

Feedback
For any feedback, questions, or suggestions, please feel free to contact at p.shamokhina@mail.ru.

