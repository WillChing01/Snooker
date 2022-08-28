# Snooker

This is a basic 2D snooker game. The game contains a functional UI as well as custom controls and the ability to choose your cue design.

To help with aiming, the game includes a shot predictor to show the trajectory of the cueball.

# Supported game modes:
- Singleplayer line up: Place the cueball anywhere and make a break from there
- LAN multiplayer: Play against another person and communicate via live text chat

# Example Gameplay (.gif) - Old Version
![clearance](https://user-images.githubusercontent.com/53403691/87889553-151ed200-ca2a-11ea-9cb6-4f8ea2755ddf.gif)

# Physics

The general motion of the balls is solved analytically rather than using a discrete timestep for greater accuracy and efficiency [1].

However the motion of a ball on the rim of a pocket is solved numerically [2].

To solve the quartic equations required for the simulation, a self-coded implementation of Algorithm 954 was used [3].

Collisions between balls were solved simulataneously using a matrix method [4]. While collisions between discrete pairs of balls could be solved using a simpler method, the matrix method was useful for providing realistic break-off shots at the start of a frame.

# References

[1] W. Leckie and M. Greenspan, "An Event-Based Pool Physics Simulator", Lecture Notes in Computer Science, pp. 247-262, 2006. Available: 10.1007/11922155_19.

[2] M. Bacon, "How balls roll off tables", American Journal of Physics, vol. 73, no. 8, pp. 722-724, 2005. Available: 10.1119/1.1947198.

[3] N. Flocke, "Algorithm 954", ACM Transactions on Mathematical Software, vol. 41, no. 4, pp. 1-24, 2015. Available: 10.1145/2699468.

[4] M. balls, A. Zita, A. Zita, J. Alexiou and s. gerbil, "Multiple colliding balls", Physics Stack Exchange, 2022. [Online]. Available: https://physics.stackexchange.com/questions/296767/multiple-colliding-balls. [Accessed: 28- Aug- 2022].
