# CarND Model Predictive Control Project

* In this project a set of way points are given for a predefined route from a simulator. Route is nothing but set of sequences in car's reference frame, where x is pointing forward and y is pointing towards left. Car's state is defined by it's current location, velocity and orientation.

* Expected output is nothing but the acceleration and turning rate from Model Predictive Controller.

* These value update the Car's state as per the given Newtonian formula in lessons. For recap, car's state is co-ordinates in global frame, orientation and velocity. Please note that the car's co-ordinates will be zero in local refrence. It's orienation will also be zero in local reference. Only velocity remains as such in the dircetion of x-local axis.

* Lf is used as the contribution linear factor of speed to the turning rate.

* Let's dt is 0.05 and N is 10, then time window is 0.7 seconds. I could find the N value from hit and trial. Larger value for N was causing drunken driving in the car. Latency of 0.05 is chosen with hit an dtrial method.

* I tried various value for N and delta delay. e.g; Setting N value of 25 cause too much wobble in the car.
